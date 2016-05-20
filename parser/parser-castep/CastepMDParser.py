import setup_paths
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from CastepCommon import get_metaInfo
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
    # def onClose_section_run(self, backend, gIndex, section):    
        # for i in range (len(temp)):
        # backend.openSection('section_single_configuration_calculation')
        # backend.addArrayValues('atom_forces', np.asarray(self.md_forces))
        # backend.closeSection('section_single_configuration_calculation',gIndex+i) 
    
    def onClose_section_md(self, backend, gIndex, section):
        temp = section ['castep_md_temperature']
        vet_veloc = section ['castep_md_cell_vectors_vel']
        velocities = section ['castep_md_veloc']
        vet = section ['castep_md_cell_vectors']
        forces = section ['castep_md_forces']
        atoms_lab = section ['castep_md_lab']
        position = section ['castep_md_positions']
        press = section ['castep_md_pressure']
        energies = section['castep_md_energies']
        stress_tens = section ['castep_md_stress_tensor']
        # frame_time =section['castep_md_time']

        
        for i in range(len(temp)):
            self.frame_temperature.append(temp[i])
  
        
        for i in range(len(press)):
            self.frame_pressure.append(press[i])
        
        for i in range(len(temp)):
            energies[i] = energies[i].split()
            energies[i] = [float(j) for j in energies[i]]
            energy_list = energies[i]

            self.total_energy.append(energy_list[0])
            
            self.hamiltonian.append(energy_list[1])
            self.kinetic.append(energy_list[2]) 
   
        if vet:
            self.cell =[]
            for i in range(len(vet)):
                vet[i] = vet[i].split()  
                vet[i] = [float(j) for j in vet[i]]                
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
            
        
        
        if stress_tens is not None:
            self.stress_tensor_value =[]
            for i in range(len(stress_tens)):
                stress_tens[i] = stress_tens[i].split()
                stress_tens[i] = [float(j) for j in stress_tens[i]]
                stress_tens_int = stress_tens[i]
                stress_tens_int = [x / 10e9 for x in stress_tens_int] #converting GPa in Pa.
                self.stress_tensor_value.append(stress_tens_int)
            self.frame_stress_tensor.append(self.stress_tensor_value)
     
        if position:
            self.at_nr = len(position)
            self.atom_position=[]
            for i in range(0, self.at_nr):
                position[i] = position[i].split()
                position[i] = [float(j) for j in position[i]]
                self.atom_position.append(position[i])
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
            startReStr =r"\s*(?P<castep_md_energies>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)\s*\<\-\-\sE\s*",
            endReStr ="/n",
            sections = ['section_md'],
            repeats = True,
            # weak = True,
            subMatchers = [
                SM (r"\s*(?P<castep_md_temperature>[-+0-9.eEdD]+)\s*\<\-\-\sT\s*"),
                SM (r"\s*(?P<castep_md_pressure>[-+0-9.eEdD]+)\s*\<\-\-\sP\s*"),
                SM (r"\s*(?P<castep_md_cell_vectors>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)\s*\<\-\-\sh\s",repeats = True),
                SM (r"\s*(?P<castep_md_cell_vectors_vel>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)\s*\<\-\-\sh[a-z]\s*",repeats = True),
                SM (r"\s*(?P<castep_md_stress_tensor>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)\s*\<\-\-\sS\s",repeats = True),
                SM(r"\s(?P<castep_md_lab>[A-Za-z]+\s*[0-9.]+)\s*(?P<castep_md_positions>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)\s*\<\-\-\sR\s*",repeats = True),
                SM(r"\s[A-Za-z]+\s*[0-9.]+\s*(?P<castep_md_veloc>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)\s*\<\-\-\sV\s*",repeats = True),
                SM(r"\s[A-Za-z]+\s*[0-9.]+\s*(?P<castep_md_forces>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)\s*\<\-\-\sF\s*",repeats = True,endReStr ="/n"),
                # SM (r"\s*(?P<castep_md_time>[+0-9.eEdD]+)\s*"),                  # ),
            

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
    cachingLevelForMetaName = { 'castep_md_energies':CachingLevel.Cache,
                                'castep_md_forces':CachingLevel.Cache,
                                'castep_md_veloc':CachingLevel.Cache,       
                                'castep_md_lab':CachingLevel.Cache,
                                'castep_md_positions':CachingLevel.Cache,
                                'castep_md_temperature':CachingLevel.Cache,
                                'castep_md_pressure':CachingLevel.Cache,
                                'castep_md_cell_vectors_vel':CachingLevel.Cache,
                                'castep_md_cell_vectors':CachingLevel.Cache,
                                'castep_md_stress_tensor':CachingLevel.Cache,
                                
                              
                                'section_run': CachingLvl,
                                'section_md': CachingLvl,
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







