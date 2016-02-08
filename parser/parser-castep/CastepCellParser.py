import setup_paths
import numpy as np
import math
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
import re, os, sys, json, logging


class CastepCellParserContext(object):

    def __init__(self, main_file_context):
        self.kp = []

    def startedParsing(self, path, parser):
        self.parser = parser

    def onClose_section_k_band(self, backend, gIndex, section):
        k_p = section['castep_store_k_path']

        k_count = len(k_p)
        k_crd = []
        for i in range(0, k_count):
            k_p[i] = k_p[i].split()
            k_p[i] = [float(j) for j in k_p[i]]
            k_p_int = k_p[i]
            k_crd.append(k_p_int)


        k_sgt_start_end = []
        for i in range(k_count-1):
            k_sgt = [ k_crd[i], k_crd[i+1] ]
            k_sgt_start_end.append(k_sgt)


        backend.addArrayValues('castep_k_path', np.asarray(k_crd))
        backend.addArrayValues('band_segm_start_end', np.asarray(k_sgt_start_end))






################################################################################
######################  MAIN PARSER STARTS HERE  ###############################
################################################################################


# main Parser
cellFileDescription = SM(name = 'root',
                         weak = True,
                         startReStr = "\! Specify a path through the Brillouin Zone to compute the band structure\.\s*",
                         subMatchers = [

                             SM(startReStr = r"\! Specify a path through the Brillouin Zone to compute the band structure\.\s*",
                                                 forwardMatch = True,
                                                 sections = ["section_k_band"],
                                                 subMatchers = [

                                                     SM(r"(?P<castep_store_k_path>[\d\.]+\s+[\d\.]+\s+[\d\.]+)", repeats = True),

                                                     ]),

])





# loading metadata from nomad-meta-info/meta_info/nomad_meta_info/fhi_aims.nomadmetainfo.json
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)
# dictionary for functions to call on close of a section metaname
onClose = {}
# parser info
parserInfo = {'name':'castep-parser', 'version': '1.0'}
# adjust caching of metadata
cachingLevelForMetaName = {
                            'castep_k_path': CachingLevel.Forward,
                            'castep_store_k_path': CachingLevel.Cache,
                           }

if __name__ == "__main__":
    mainFunction(cellFileDescription, metaInfoEnv, parserInfo, superContext = CastepCellParserContext(), cachingLevelForMetaName = cachingLevelForMetaName, onClose = onClose,
                 defaultSectionCachingLevel = True)
