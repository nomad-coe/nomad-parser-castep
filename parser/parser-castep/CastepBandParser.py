import setup_paths
import numpy as np
import math
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
import re, os, sys, json, logging


class CastepParserContext(object):

    def __init__(self):
        self.k_count                    = 0
        self.k_nr                       = 0
        self.e_nr                       = 0
        self.eigenvalues_kpoints        = []
        self.eigenvalues_eigenvalues    = []



# Storing the k point coordinates
    def onClose_castep_section_k_points(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""

# Processing k points (given in fractional coordinates)
        #get cached values of castep_store_k_points
        k_st = section.simpleValues['castep_store_k_points']
        self.k_count = len(k_st)
        self.k_nr   += 1
        for i in range(0, self.k_count):
            k_st[i] = k_st[i].split()
            k_st[i] = [float(j) for j in k_st[i]]
            k_st_int = k_st[i]
            self.eigenvalues_kpoints.append(k_st_int)


# Storing the eigenvalues
    def onClose_castep_section_eigenvalues(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""

        #get cached values of castep_store_k_points
        e_st = section.simpleValues['castep_store_eigenvalues']
        self.e_nr = len(e_st)
        self.eigenvalues_eigenvalues.append(e_st)



# Keeping only the last arrays update
    def onClose_section_eigenvalues(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
        backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.eigenvalues_kpoints))
        backend.addArrayValues('eigenvalues_eigenvalues', np.asarray(self.eigenvalues_eigenvalues))
# Backend add the number of k points and eigenvalues
        backend.addValue('eigenvalues_kpoints_number', self.k_nr)
        backend.addValue('eigenvalues_eigenvalues_number', self.e_nr)








################################################################################
######################  MAIN PARSER STARTS HERE  ###############################
################################################################################

# main Parser
mainFileDescription = SM(name = 'root',
                         weak = True,
                         startReStr = "",
                         subMatchers = [

                            SM(name = "newRun",
                            startReStr = r"",
                            #repeats = True,
                            required = True,
                            forwardMatch = True,
                            sections   = ["section_run"],
                            subMatchers = [


                                SM(startReStr = "Number of k-points\s*",
                                   forwardMatch = True,
                                   sections = ["section_system_description"],
                                   subMatchers = [

                                       SM(startReStr = "Number of k-points\s*",
                                          forwardMatch = True,
                                          sections = ["section_single_configuration_calculation"],
                                          subMatchers = [

                                              SM(startReStr = "Number of k-points\s*",
                                                 forwardMatch = True,
                                                 sections = ["section_eigenvalues_group"],
                                                 subMatchers = [

                                                     SM(startReStr = "Number of k-points\s*",
                                                        sections = ["section_eigenvalues"],
                                                        forwardMatch = True,
                                                        subMatchers = [

                                                         SM(startReStr = "K\-point\s*[0-9]+\s*",
                                                            sections = ["castep_section_k_points"],
                                                            forwardMatch = True,
                                                            repeats = True,
                                                            subMatchers = [

                                                                SM(r"K\-point\s*[0-9]+\s*(?P<castep_store_k_points> [-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                                                   repeats = True),

                                                                    SM(name = 'Eigen',
                                                                       startReStr = r"Spin component\s*[0-9]+\s*",
                                                                       sections = ['castep_section_eigenvalues'],
                                                                       repeats = True,
                                                                       subMatchers = [

                                                                          SM(r"\s*(?P<castep_store_eigenvalues> [-\d\.]+)",
                                                                             repeats = True)

                                                                       ]), # CLOSING castep_section_eigenvalues


                                                            ]), # CLOSING castep_section_k_points


                                                        ]), # CLOSING section_eigenvalues



                                                ]), # CLOSING section_eigenvalues_group


                                          ]), # CLOSING section_single_configuration_calculation


                                ]), # CLOSING section_system_description


    ]), # CLOSING SM newRun


])






# loading metadata from nomad-meta-info/meta_info/nomad_meta_info/fhi_aims.nomadmetainfo.json
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)
# dictionary for functions to call on close of a section metaname
onClose = {}
# parser info
parserInfo = {'name':'castep-parser', 'version': '1.0'}
# adjust caching of metadata
cachingLevelForMetaName = {'eigenvalues_kpoints': CachingLevel.Forward,
                           'eigenvalues_eigenvalues': CachingLevel.Forward,
                           'eigenvalues_kpoints_number': CachingLevel.Forward,
                           'eigenvalues_eigenvalues_number': CachingLevel.Forward,
                           'castep_store_k_points': CachingLevel.Cache,
                           'castep_store_eigenvalues': CachingLevel.Cache,
                           }

if __name__ == "__main__":
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo, superContext = CastepParserContext(), cachingLevelForMetaName = cachingLevelForMetaName, onClose = onClose)
