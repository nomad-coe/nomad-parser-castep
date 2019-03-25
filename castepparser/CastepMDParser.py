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



############################################################
# This is the parser for the *.md file of CASTEP.
############################################################

logger = logging.getLogger("nomad.CastepMDParser")

class CastepMDParserContext(object):
    """Context for parsing CASTEP *.md file.


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
        self.frame_temperature =[]
        self.frame_pressure = []
        self.total_energy = []
        self.hamiltonian = []
        self.kinetic = []
        self.md_forces =[]
        self.md_veloc =[]
        self.frame_atom_label =[]
        self.total_velocities =[]
        self.total_forces =[]
        self.atom_label = []
        self.frame_stress_tensor =[]
        self.total_positions =[]
        self.frame_cell=[]
        self.vector_velocities=[]
        self.frame_time =[]
    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts and the compiled parser is obtained.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser

    def onClose_x_castep_section_md(self, backend, gIndex, section):
        temp = section ['x_castep_md_temperature']
        vet_veloc = section ['x_castep_md_cell_vectors_vel']
        velocities = section ['x_castep_md_veloc']
        vet = section ['x_castep_md_cell_vectors']
        forces = section ['x_castep_md_forces']
        atoms_lab = section ['x_castep_md_lab']
        position = section ['x_castep_md_positions']
        press = section ['x_castep_md_pressure']
        energies = section['x_castep_md_energies']
        stress_tensor = section ['x_castep_md_stress_tensor']

        Hr_J_converter = float(4.35974e-18)
        HrK_to_K_coverter= float(3.1668114e-6)
        evAtoN = float(1.6021766e-9)
        a_to_m = 1e-10



        for i in range(len(temp)):
            temp[i] = temp[i]/HrK_to_K_coverter
            self.frame_temperature.append(temp[i])


        for i in range(len(press)):
            press[i] = press[i] / 10e9
            self.frame_pressure.append(press[i])

        for i in range(len(temp)):
            energies[i] = energies[i].split()
            energies[i] = [float(j) for j in energies[i]]
            energy_list = energies[i]
            energy_list = [x * Hr_J_converter for x in energy_list]
            self.total_energy.append(energy_list[0])
            self.hamiltonian.append(energy_list[1])
            self.kinetic.append(energy_list[2])

        if vet:
            self.cell =[]
            for i in range(len(vet)):
                vet[i] = vet[i].split()
                vet[i] = [float(j) for j in vet[i]]
                vet[i] = [ii * a_to_m  for ii in vet[i]]
                vet_list = vet[i]

                self.cell.append(vet_list)
        self.frame_cell.append(self.cell)


        if vet_veloc:
            self.vet_vel =[]
            for i in range(len(vet_veloc)):
                vet_veloc[i] = vet_veloc[i].split()
                vet_veloc[i] = [float(k) for k in vet_veloc[i]]

                v_vet = vet_veloc[i]
                self.vet_vel.append(v_vet)
            self.vector_velocities.append(self.vet_vel)

        if stress_tensor is not None:
            self.stress_tensor_value =[]
            for s in stress_tensor:
                s = s.split()
                s = [float(k) for k in s]
                stress_tens_int = s
                stress_tens_int = [x / 10e9 for x in stress_tens_int]
                self.stress_tensor_value.append(stress_tens_int)
            self.frame_stress_tensor.append(self.stress_tensor_value)

        if position:
            self.at_nr = len(position)
            self.atom_position=[]
            for i in range(0, self.at_nr):
                position[i] = position[i].split()
                position[i] = [float(j) for j in position[i]]
                position[i] = [ii * a_to_m  for ii in position[i]]
                pos_list = position[i]
                self.atom_position.append(pos_list)
            self.total_positions.append(self.atom_position)


        if velocities is not None:
            self.md_veloc =[]
            for j in range(len(velocities)):
                velocities[j] = velocities[j].split()
                velocities[j] = [float(k) for k in velocities[j]]
                v_st_int = velocities[j]
                self.md_veloc.append(v_st_int)
            self.total_velocities.append(self.md_veloc)

        if forces is not None:

            self.md_forces = []
            for f in forces:
                f = f.split()
                f = [float(k) for k in f]
                f_st_int = f
                f_st_int = [x * evAtoN for x in f_st_int]
                self.md_forces.append(f_st_int)
            self.total_forces.append(self.md_forces)







def build_CastepMDFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the *.md file of CASTEP.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses *.md file of CASTEP.
    """
    return SM (name = 'Root1',
        startReStr = "",
        sections = ['section_run'],
        forwardMatch = True,
        weak = True,
        subMatchers = [
            SM (name = 'Root2',
            startReStr =r"\s*(?P<x_castep_md_energies>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sE\s*",
            endReStr ="/n",
            sections = ['x_castep_section_md'],
            repeats = True,
            subMatchers = [
                SM (r"\s*(?P<x_castep_md_temperature>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sT\s*"),
                SM (r"\s*(?P<x_castep_md_pressure>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sP\s*"),
                SM (r"\s*(?P<x_castep_md_cell_vectors>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sh\s",repeats = True),
                SM (r"\s*(?P<x_castep_md_cell_vectors_vel>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sh[a-z]\s*",repeats = True),
                SM (r"\s*(?P<x_castep_md_stress_tensor>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sS\s",repeats = True),
                SM(r"\s(?P<x_castep_md_lab>[A-Za-z]+\s*[0-9.]+)\s*(?P<x_castep_md_positions>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sR\s*",repeats = True),
                SM(r"\s[A-Za-z]+\s*[0-9.]+\s*(?P<x_castep_md_veloc>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sV\s*",repeats = True),
                SM(r"\s[A-Za-z]+\s*[0-9.]+\s*(?P<x_castep_md_forces>[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s*[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)\s*\<\-\-\sF\s*",repeats = True,endReStr ="/n"),

            ]),
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
    cachingLevelForMetaName = { 'x_castep_md_energies':CachingLevel.Cache,
                                'x_castep_md_forces':CachingLevel.Cache,
                                'x_castep_md_veloc':CachingLevel.Cache,
                                'x_castep_md_lab':CachingLevel.Cache,
                                'x_castep_md_positions':CachingLevel.Cache,
                                'x_castep_md_temperature':CachingLevel.Cache,
                                'x_castep_md_pressure':CachingLevel.Cache,
                                'x_castep_md_cell_vectors_vel':CachingLevel.Cache,
                                'x_castep_md_cell_vectors':CachingLevel.Cache,
                                'x_castep_md_stress_tensor':CachingLevel.Cache,


                                'section_run': CachingLvl,
                                'x_castep_section_md': CachingLvl,
                               # 'section_single_configuration_calculation': CachingLvl,
                              }
    # Set all band metadata to Cache as they need post-processsing.
    # for name in metaInfoEnv.infoKinds:
    #     if name.startswith('castep_'):
    #         cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main(CachingLvl):
    """Main function.

    Set up everything for the parsing of the CASTEP *.cell file and run the parsing.

    Args:
        CachingLvl: Sets the CachingLevel for the sections k_band, run, and single_configuration_calculation.
            This allows to run the parser without opening new sections.
    """
    # get band.out file description
    CastepMDFileSimpleMatcher = build_CastepMDFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name':'castep-md-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv, CachingLvl)
    # start parsing
    mainFunction(mainFileDescription = CastepMDFileSimpleMatcher,
                 metaInfoEnv = metaInfoEnv,
                 parserInfo = parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = CastepMDParserContext())

if __name__ == "__main__":
    main(CachingLevel.Forward)







