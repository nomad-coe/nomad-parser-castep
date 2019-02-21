# Copyright 2015-2018 Martina Stella, Massimo Riello, Fawzi Mohamed, Ankit Kariryaa
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from builtins import range
from builtins import object
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from castepparser.CastepCommon import get_metaInfo
import logging, os, re, sys



############################################################
# This is the parser for the *.cell file of CASTEP.
############################################################

logger = logging.getLogger("nomad.CastepCellParser")

class CastepCellParserContext(object):
    """Context for parsing CASTEP *.cell file.


    The onClose_ functions allow processing and writing of cached values after a section is closed.
    They take the following arguments:
        backend: Class that takes care of wrting and caching of metadata.
        gIndex: Index of the section that is closed.
        section: The cached values and sections that were found in the section that is closed.
    """
    def __init__(self, writeMetaData = True):
        """Args:
            writeMetaData: Deteremines if metadata is written or stored in class attributes.
        """
        self.writeMetaData = writeMetaData

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts and the compiled parser is obtained.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        # get unit from metadata for band energies
        # allows to reset values if the same superContext is used to parse different files
        self.band_energies = None
        self.band_k_points = None
        self.band_occupations = None
        self.klabel_seg = []
        self.k_crd = []
        self.k_sgt_start_end = []

    def onClose_section_k_band(self, backend, gIndex, section):
        """Trigger called when section_k_band is closed.

        Store the parsed values and write them if writeMetaData is True.
        """
        k_p = section['x_castep_store_k_path']

        k_count = len(k_p)
        self.k_crd = []
        for i in range(0, k_count):
            k_p[i] = k_p[i].split()
            k_p[i] = [float(j) for j in k_p[i]]
            k_p_int = k_p[i]
            self.k_crd.append(k_p_int)


        self.k_sgt_start_end = []
        for i in range(k_count-1):
            k_sgt = [ self.k_crd[i], self.k_crd[i+1] ]
            self.k_sgt_start_end.append(k_sgt)


        backend.addArrayValues('x_castep_k_path', np.asarray(self.k_crd))
        self.klabel_seg = section['x_castep_store_k_label']

        # backend.addArrayValues('band_segm_start_end', np.asarray(self.k_sgt_start_end))
        # backend.addValue('number_of_k_point_segments', len(self.k_sgt_start_end))


def build_CastepCellFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the *.cell file of CASTEP.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses *.cell file of CASTEP.
    """
    return SM (name = 'Root1',
        startReStr = "",
        sections = ['section_run'],
        forwardMatch = True,
        weak = True,
        subMatchers = [
        SM (name = 'Root2',
            startReStr = "",
            sections = ['section_single_configuration_calculation'],
            forwardMatch = True,
            weak = True,
            subMatchers = [

            SM(startReStr = r"\s*\%block bs\_kpoint\_path\s*",
                  sections = ['section_k_band'],
                  forwardMatch = True,
                  subMatchers = [

                     SM (r"(?P<x_castep_store_k_path>[\d\.]+\s+[\d\.]+\s+[\d\.]+)\s*\!\s(?P<x_castep_store_k_label>[A-Z])", repeats = True)
                                 ]),

            SM(startReStr = r"\s*\%BLOCK BS\_KPOINT\_PATH\s*",
                  sections = ['section_k_band'],
                  forwardMatch = True,
                  subMatchers = [

                      SM (r"(?P<x_castep_store_k_path>[\d\.]+\s+[\d\.]+\s+[\d\.]+)\s*\!\s(?P<x_castep_store_k_label>[A-Z])", repeats = True)

                                 ]),


            # SM(startReStr = r"\s*\%block kpoint\_list\s*",
            #       sections = ['section_k_band'],
            #       forwardMatch = True,
            #       subMatchers = [

            #          SM (r"\s*(?P<x_castep_store_k_path>[\d\.]+\s+[\d\.]+\s+[\d\.]+[\d\.]+)", repeats = True)
            #                      ]),

            # SM(startReStr = r"\s*\%BLOCK KPOINT\_LIST\s*",
            #       sections = ['section_k_band'],
            #       forwardMatch = True,
            #       subMatchers = [

            #           SM (r"\s*(?P<x_castep_store_k_path>[\d\.]+\s+[\d\.]+\s+[\d\.]+[\d\.]+)", repeats = True)

            #                      ]),

            SM(startReStr = r"\%block kpoint\_list\s*",
                  sections = ['section_k_band'],
                  forwardMatch = True,
                  subMatchers = [

                     SM (r"\s*(?P<x_castep_store_k_path>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)\s*\!\s(?P<x_castep_store_k_label>[A-Z])", repeats = True)
                                 ]),

            # SM(startReStr = r"\%BLOCK KPOINT\_LIST\s*",
            #       sections = ['section_k_band'],
            #       forwardMatch = True,
            #       subMatchers = [

            #           SM (r"\s*(?P<x_castep_store_k_path>[\d\.]+\s+[\d\.]+\s+[\d\.]+[\d\.]+)", repeats = True)

            #                      ]),
            #SM (name = 'Root3',
            #    startReStr = r"\s*\%block bs\_kpoint\_path\s*",
            #    sections = ['section_k_band'],
            #    forwardMatch = True,
            #    weak = True,
            #    subMatchers = [
            #    SM (r"(?P<castep_store_k_path>[\d\.]+\s+[\d\.]+\s+[\d\.]+)", repeats = True)
            #    ]),

            ])

        ])


def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
        CachingLvl: Sets the CachingLevel for the sections k_band, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
                               'section_k_band': CachingLvl,
                               'section_run': CachingLvl,
                               'section_single_configuration_calculation': CachingLvl,
                              }
    # Set all band metadata to Cache as they need post-processsing.
    for name in metaInfoEnv.infoKinds:
        if name.startswith('x_castep_'):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main(CachingLvl):
    """Main function.

    Set up everything for the parsing of the CASTEP *.cell file and run the parsing.

    Args:
        CachingLvl: Sets the CachingLevel for the sections k_band, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.
    """
    # get band.out file description
    CasteCellFileSimpleMatcher = build_CastepCellFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name':'castep-cell-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv, CachingLvl)
    # start parsing
    mainFunction(mainFileDescription = CasteCellFileSimpleMatcher,
                 metaInfoEnv = metaInfoEnv,
                 parserInfo = parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = CastepCellParserContext())

if __name__ == "__main__":
    main(CachingLevel.Forward)







