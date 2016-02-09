import setup_paths
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from CastepCommon import get_metaInfo
import CastepCellParser_1
import logging, os, re, sys





################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
######################  PARSER CONTEXT CLASS  ##################################################################################################################
################################################################################################################################################################
############################################################################ CASTEP.Parser Version 1.0 #########################################################
################################################################################################################################################################

class CastepParserContext(object):

    def __init__(self):
        self.cell                              = []


    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        # start with -1 since zeroth iteration is the initialization
        self.scfIterNr = -1


    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Get compiled parser, filename and metadata.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        self.fName = fInName
        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()


    def onClose_section_k_band(self, backend, gIndex, section):
        """Trigger called when section_k_band is closed.

        Band structure is parsed from external band.out files.
        """

        cellSuperContext = CastepCellParser_1.CastepCellParserContext(False)
        cellParser = AncillaryParser(
            fileDescription = CastepCellParser_1.build_CastepCellFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = CastepCellParser_1.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = cellSuperContext)

        bFile = "Si2.cell"
        dirName = os.path.dirname(os.path.abspath(self.fName))
        fName = os.path.normpath(os.path.join(dirName, bFile))

        with open(fName) as fIn:
            cellParser.parseFile(fIn)

        k_start_end = cellSuperContext.k_sgt_start_end








def build_CastepMainFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the main file of FHI-aims.

    First, several subMatchers are defined, which are then used to piece together
    the final SimpleMatcher.
    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses main file of FHI-aims.
    """


    ########################################
    # submatcher for band structure
    bandStructureSubMatcher = SM (name = 'BandStructure',
        startReStr = r"\s*\+\s*B A N D   S T R U C T U R E   C A L C U L A T I O N\s*",
        sections = ['section_k_band'],
        subMatchers = [


            # First spin channel
            SM(startReStr = "\s*\+\s*Spin\=1\skpt\=\s*",
               sections = ["castep_section_k_band"],
               forwardMatch = True,
               subMatchers = [

                  SM(startReStr = "\s*\+\s*Spin\=1\skpt\=\s*",
                     sections = ["castep_section_k_points"],
                     forwardMatch = True,
                     repeats = True,
                     subMatchers = [
                        # Matching k points
                        SM(r"\s*\+\s*Spin\=1\s*kpt\=\s*[0-9]+\s*\((?P<castep_store_k_points>\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)\)\s*",
                           repeats = True),

                           SM(name = 'Eigen_1',
                              startReStr = r"\s*\+\s*\+\s*",
                              sections = ['castep_section_eigenvalues'],
                              repeats = True,
                              subMatchers = [
                                 # Matching eigenvalues
                                 SM(r"\s*\+\s*[0-9]+\s*(?P<castep_store_eigenvalues>\s+[-\d\.]+)",
                                    repeats = True)

                                             ]), # CLOSING castep_section_eigenvalues


                                    ]), # CLOSING castep_section_k_points


                              ]), # CLOSING 1st section_eigenvalues

            # Second spin channel
            SM(startReStr = "\s*\+\s*Spin\=2\skpt\=\s*",
               sections = ["castep_section_k_band"],
               forwardMatch = True,
               subMatchers = [

                  SM(startReStr = "\s*\+\s*Spin\=2\skpt\=\s*",
                     sections = ["castep_section_k_points_1"],
                     forwardMatch = True,
                     repeats = True,
                     subMatchers = [
                        # Matching k points
                        SM(r"\s*\+\s*Spin\=2\s*kpt\=\s*[0-9]+\s*\((?P<castep_store_k_points_1>\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)\)\s*",
                           repeats = True),

                           SM(name = 'Eigen_2',
                              startReStr = r"\s*\+\s*\+\s*",
                              sections = ['castep_section_eigenvalues_1'],
                              repeats = True,
                              subMatchers = [
                                 # Matching eigenvalues
                                 SM(r"\s*\+\s*[0-9]+\s*(?P<castep_store_eigenvalues_1>\s+[-\d\.]+)",
                                    repeats = True)

                                             ]), # CLOSING castep_section_eigenvalues_1


                                    ]), # CLOSING castep_section_k_points_1


                              ]), # CLOSING 2nd section_eigenvalues



        ])


    ########################################
    # return main Parser
    return SM (name = 'Root',
        startReStr = "",
        forwardMatch = True,
        weak = True,
        subMatchers = [
        SM (name = 'NewRun',
            startReStr = r"\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*",
            required = True,
            forwardMatch = True,
            sections = ['section_run'],
            subMatchers = [

               SM(startReStr = r"\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*",
                  forwardMatch = True,
                  sections = ["section_single_configuration_calculation"],
                  subMatchers = [

                                 bandStructureSubMatcher

                                 ])

                           ])

        ])







def get_cachingLevelForMetaName(metaInfoEnv):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = { }

    # Set all controlIn and controlInOut metadata to Cache to capture multiple occurrences of keywords and
    # their last value is then written by the onClose routine in the FhiAimsParserContext.
    # Set all geometry metadata to Cache as all of them need post-processsing.
    # Set all eigenvalue related metadata to Cache.
    for name in metaInfoEnv.infoKinds:
        if (   name.startswith('castep_store_')
            or name.startswith('fhi_aims_geometry_')
            or name.startswith('fhi_aims_eigenvalue_')
            or name.startswith('fhi_aims_section_eigenvalues_')
           ):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName

def main():
    """Main function.

    Set up everything for the parsing of the FHI-aims main file and run the parsing.
    """
    # get main file description
    CastepMainFileSimpleMatcher = build_CastepMainFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/fhi_aims.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name':'castep-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv)
    # start parsing
    mainFunction(mainFileDescription = CastepMainFileSimpleMatcher,
                 metaInfoEnv = metaInfoEnv,
                 parserInfo = parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = CastepParserContext(),
                 defaultSectionCachingLevel = False)

if __name__ == "__main__":
    main()








