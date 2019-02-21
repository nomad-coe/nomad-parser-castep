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

from __future__ import division
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




###############################################################
# This is the parser for the *.band file of CASTEP.


# NB: this parser store self consistent eigenvalues and
#     relative k points for single configuration calculations

###############################################################

logger = logging.getLogger("nomad.CastepBandParser")

class CastepBandParserContext(object):
    """Context for parsing CASTEP *.band file.


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
        self.k_count                    = 0
        self.k_nr                       = 0
        self.e_nr                       = 0
        self.eigenvalues_kpoints        = []
        self.eigenvalues_values         = []

        self.e_spin_1                   = []
        self.e_spin_2                   = []
        self.n_spin                     = 0



# Reading the number of spins
    def onClose_x_castep_section_spin_number(self, backend, gIndex, section):

        self.n_spin = section['x_castep_spin_number']


# Storing the k point coordinates
    def onClose_x_castep_section_scf_k_points(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""

# Processing k points (given in fractional coordinates)
        #get cached values of x_castep_store_k_points
        k_st = section['x_castep_store_scf_k_points']
        self.k_count = len(k_st)
        self.k_nr   += 1
        for i in range(0, self.k_count):
            k_st[i] = k_st[i].split()
            k_st[i] = [float(j) for j in k_st[i]]
            k_st_int = k_st[i]
            self.eigenvalues_kpoints.append(k_st_int)


# Storing the eigenvalues
    def onClose_x_castep_section_scf_eigenvalues(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
        Ha_to_J = 4.35974e-18
        #get cached values of castep_store_k_points
        e_st = section['x_castep_store_scf_eigenvalues']

        e_st_0 = e_st
        e_st_0 = [x * Ha_to_J for x in e_st_0]

        def split_list(lista):
            half = len(lista) // 2
            return lista[:half], lista[half:]

        e_st_1, e_st_2 = split_list(e_st)
        e_st_1 = [x * Ha_to_J for x in e_st_1]
        e_st_2 = [x * Ha_to_J for x in e_st_2]

        if self.n_spin[0] == 1:
            self.e_nr = len(e_st_0)
            self.e_spin_1.append(e_st_0)
            self.e_spin_2 = []

        else:
            self.e_nr = len(e_st_1)
            self.e_spin_1.append(e_st_1)
            self.e_spin_2.append(e_st_2)







################################################################################
######################  MAIN PARSER STARTS HERE  ###############################
################################################################################


def build_CastepBandFileSimpleMatcher():
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

               SM(startReStr = r"Number\sof\sk\-points\s*[0-9]+\s*",
                  sections = ['x_castep_section_spin_number'],
                  forwardMatch = True,
                  subMatchers = [

                     SM(r"Number\sof\sspin\scomponents\s*(?P<x_castep_spin_number>[1-2]+)")

                                 ]),

               SM(startReStr = r"K\-point\s*[0-9]+\s*",
                  sections = ["x_castep_section_scf_k_points"],
                  forwardMatch = True,
                  repeats = True,
                  subMatchers = [

                     SM(r"K\-point\s*[0-9]+\s*(?P<x_castep_store_scf_k_points> [-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                        repeats = True),

                     SM(name = 'Eigen',
                        startReStr = r"Spin component\s*1\s*",
                        #endReStr = r"Spin component\s*2\s*",
                        sections = ['x_castep_section_scf_eigenvalues'],
                        repeats = True,
                        subMatchers = [

                           SM(r"\s*(?P<x_castep_store_scf_eigenvalues> [-\d\.]+)",
                              repeats = True)

                                       ]), # CLOSING castep_section_scf_eigenvalues


                                 ]) # CLOSING castep_section_k_points


    ]), # CLOSING section_single_configuration_calculation
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
    CasteBandFileSimpleMatcher = build_CastepBandFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name':'castep-cell-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv, CachingLvl)
    # start parsing
    mainFunction(mainFileDescription = CasteBandFileSimpleMatcher,
                 metaInfoEnv = metaInfoEnv,
                 parserInfo = parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = CastepBandParserContext())

if __name__ == "__main__":
    main(CachingLevel.Forward)


