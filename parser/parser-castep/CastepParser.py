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
        """

        self.k_nr                              = 0
        self.e_nr                              = 0
        self.k_count_1                         = 0
        self.k_nr_1                            = 0
        self.e_nr_1                            = 0
        self.castep_band_kpoints               = []
        self.castep_band_energies              = []
        self.castep_band_kpoints_1             = []
        self.castep_band_energies_1            = []
        self.k_path_nr                         = 0

        self.band_en = []


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


######################################################################################
###################### Storing k points and band energies ############################
############################# FIRST SPIN CHANNEL #####################################
######################################################################################

# Storing the k point coordinates (SPIN 1)
    def onClose_castep_section_k_points(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
# Processing k points (given in fractional coordinates)
        #get cached values of castep_store_k_points
        k_st = section['castep_store_k_points']
        self.k_count = len(k_st)
        self.k_nr   += 1
        for i in range(0, self.k_count):
            k_st[i] = k_st[i].split()
            k_st[i] = [float(j) for j in k_st[i]]
            k_st_int = k_st[i]
            self.castep_band_kpoints.append(k_st_int)



# Storing the eigenvalues (SPIN 1)
    def onClose_castep_section_eigenvalues(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
        #get cached values of castep_store_k_points
        e_st = section['castep_store_eigenvalues']
        self.e_nr = len(e_st)
        self.castep_band_energies.append(e_st)


######################################################################################
###################### Storing k points and band energies ############################
############################# SECOND SPIN CHANNEL ####################################
######################################################################################

# Storing the k point coordinates (SPIN 2)
    def onClose_castep_section_k_points_1(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
# Processing k points (given in fractional coordinates)
        #get cached values of castep_store_k_points
        k_st_1 = section['castep_store_k_points_1']
        self.k_count_1 = len(k_st_1)
        self.k_nr_1   += 1
        for i in range(0, self.k_count_1):
            k_st_1[i] = k_st_1[i].split()
            k_st_1[i] = [float(j) for j in k_st_1[i]]
            k_st_1_int = k_st_1[i]
            self.castep_band_kpoints_1.append(k_st_1_int)

        self.k_nr_1 = self.k_nr  # clean double counting



# Storing the eigenvalues (SPIN 2)
    def onClose_castep_section_eigenvalues_1(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
        #get cached values of castep_store_k_points
        e_st_1 = section['castep_store_eigenvalues_1']
        self.e_nr_1 = len(e_st_1)
        self.castep_band_energies_1.append(e_st_1)

        self.e_nr_1 = self.e_nr

######################################################################################
######################################################################################
######################################################################################


    def onClose_section_k_band(self, backend, gIndex, section):
        """Trigger called when section_k_band is closed.

           The k path coordinates are extracted from the *.cell input file.
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

        self.k_start_end = cellSuperContext.k_sgt_start_end
        self.k_path_nr = len(self.k_start_end)


        if self.castep_band_energies_1 != []:
            for i in range(self.k_nr):
                a = [ self.castep_band_energies[i], self.castep_band_energies_1[i] ]
                self.band_en.append(a)
        else:
            self.band_en = self.castep_band_energies


        def get_last_index(el, check):  # function that returs end index for each k path
            found = None
            for i, next in enumerate(check):
                if next == el:
                    found = i

            assert found != None
            return found


        path_end_index = []
        for i in range(self.k_path_nr):
            #print self.k_start_end[i][1]
            boundary = self.k_start_end[i][1]
            a = get_last_index(boundary, self.castep_band_kpoints)
            path_end_index.append(a)

        path_end_index = [0] + path_end_index


        k_point_path = []
        for i in range(self.k_path_nr):
            a = self.castep_band_kpoints[ path_end_index[i] : path_end_index[i+1]+1 ]
            k_point_path.append(a)


        band_en_path = []
        for i in range(self.k_path_nr):
            a = self.band_en[ path_end_index[i] : path_end_index[i+1]+1 ]
            band_en_path.append(a)


        backend.addArrayValues('band_k_points', np.asarray(k_point_path))
        backend.addArrayValues('band_energies', np.asarray(band_en_path))




################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
######################  MAIN PARSER STARTS HERE  ###############################################################################################################
################################################################################################################################################################
############################################################################ CASTEP.Parser Version 1.0 #########################################################
################################################################################################################################################################


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
        if (   name.startswith('castep_store_')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main():
    """Main function.

    Set up everything for the parsing of the CASTEP main file and run the parsing.
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








