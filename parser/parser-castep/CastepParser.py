import setup_paths
import numpy as np
import math
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from CastepCommon import get_metaInfo
import CastepCellParser
import CastepBandParser
import logging, os, re, sys


################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
######################  PARSER CONTEXT CLASS  ##################################################################################################################
################################################################################################################################################################
############################################################################ CASTEP.Parser Version 1.0 #########################################################
################################################################################################################################################################

logger = logging.getLogger("nomad.CastepParser")

class CastepParserContext(object):

    def __init__(self):
        """ Initialise variables used within the current superContext """
        self.functionals                       = []
        self.relativistic                      = []
        self.cell                              = []
        self.at_nr                             = 0
        self.atom_type_mass                    = []
        self.atom_label                        = []
        self.atom_forces                       = []
        self.stress_tensor_value               = []
        self.castep_atom_position              = []
        self.atom_position                     = []
        self.a                                 = []
        self.b                                 = []
        self.c                                 = []
        self.alpha                             = []
        self.beta                              = []
        self.gamma                             = []
        self.volume                            = 0
        self.time_calc                         = 0
        self.time_calculation                  = []
        self.energy_total_scf_iteration_list   = []
        self.scfIterNr                         = []

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
        self.band_en                           = []

        self.e_spin_1                          = []
        self.e_spin_2                          = []


    def initialize_values(self):
        """ Initializes the values of variables in superContexts that are used to parse different files """
        self.pippo = None


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
    
    def onClose_section_topology(self, backend, gIndex, section):
    # Processing the atom masses and type
        #get cached values of castep_store_atom_mass and type
        mass = section['castep_store_atom_mass']
        name = section['castep_store_atom_name']
        for i in range(len(mass)):
            #converting mass in Kg from AMU
            mass[i] = float(mass[i]) * 1.66053904e-27

            backend.openSection('section_atom_type')
            backend.addValue('atom_type_mass', float(mass[i]))

            backend.addValue('atom_type_name', name[i])
            backend.closeSection('section_atom_type', gIndex+i)
    

# Translating the XC functional name to the NOMAD standard
    def onClose_castep_section_functionals(self, backend, gIndex, section):
        """When all the functional definitions have been gathered, matches them
        with the nomad correspondents and combines into one single string which
        is put into the backend.
        """
        # Get the list of functional and relativistic names
        functional_names = section["castep_functional_name"]
        relativistic_names = section["castep_relativity_treatment_scf"]

        # Define a mapping for the functionals
        functional_map = {
            " Perdew Burke Ernzerhof": "GGA_C_PBE_GGA_X_PBE",
            " Local Density Approximation": "LDA_C_PZ_LDA_X_PZ",
            " Perdew Wang (1991)": "GGA_C_PW91_GGA_X_PW91",
            " revised Perdew Burke Ernzerhof": "GGA_X_RPBE",
            " PBE for solids (2008)": "GGA_X_PBE_SOL",          
            " hybrid B3LYP"       :"HYB_GGA_XC_B3LYP5",
            " Hartree-Fock"       :"HF_X", 
            " Hartree-Fock + Local Density Approximation":"HF_X_LDA_C_PW",
            " hybrid HSE03"              :"HYB_GGA_XC_HSE03",
            " hybrid HSE06"              :"HYB_GGA_XC_HSE06",
            " hybrid PBE0"               :"GGA_C_PBE_GGA_X_PBE", 
            " PBE with Wu-Cohen exchange" :"GGA_C_PBE_GGA_X_WC",    
        }
       
        # Define a mapping for the relativistic treatments
        relativistic_map = {
            " Koelling-Harmon": "scalar_relativistic"
        }

        # Match each castep functional name and sort the matches into a list
        self.functionals = []

        for name in functional_names:
            match = functional_map.get(name)
            if match:
                self.functionals.append(match)
        self.functionals = "_".join(sorted(self.functionals))

        # Match each castep relativity treatment name and sort the matches into a list
        self.relativistic = []

        for name in relativistic_names:
            match = relativistic_map.get(name)
            if match:
                self.relativistic.append(match)
        self.relativistic = "_".join(sorted(self.relativistic))


# Here we add info about the XC functional and relativistic treatment
    def onClose_section_method(self, backend, gIndex, section):

        # Push the functional string into the backend
        backend.addValue('XC_functional', self.functionals)
        # Push the relativistic treatment string into the backend
        backend.addValue('relativity_method', self.relativistic)
        backend.addValue('XC_method_current', self.functionals+'_'+self.relativistic)


# Here we add basis set name and kind for the plane wave code
    def onClose_section_basis_set_cell_associated(self, backend, gIndex, section):
        ecut_str = section['castep_basis_set_plan_wave_cutoff']
        self.ecut = float(ecut_str[0])
        eVtoRy = 0.073498618
        ecut_str_name = int(round(eVtoRy*self.ecut))

        basis_set_kind = 'plane_waves'
        basis_set_name = 'PW_'+str(ecut_str_name)
        backend.addValue('basis_set_plan_wave_cutoff', self.ecut)
        backend.addValue('basis_set_cell_associated_kind', basis_set_kind)
        backend.addValue('basis_set_cell_associated_name', basis_set_name)


# Storing the unit cell
    def onClose_castep_section_cell(self, backend, gIndex, section):
        """trigger called when _castep_section_cell is closed"""
        # get cached values for castep_cell_vector
        vet = section['castep_cell_vector']

        vet[0] = vet[0].split()
        vet[0] = [float(i) for i in vet[0]]

        vet[1] = vet[1].split()
        vet[1] = [float(i) for i in vet[1]]

        vet[2] = vet[2].split()
        vet[2] = [float(i) for i in vet[2]]

        self.cell.append(vet[0])
        self.cell.append(vet[1])
        self.cell.append(vet[2]) # Reconstructing the unit cell vector by vector


# Here we recover the unit cell dimensions (both magnitudes and angles) (useful to convert fractional coordinates to cartesian)
    def onClose_castep_section_atom_position(self, backend, gIndex, section):
        """trigger called when _castep_section_atom_position is closed"""
        # get cached values for cell magnitudes and angles
        self.a = section['castep_cell_length_a']
        self.b = section['castep_cell_length_b']
        self.c = section['castep_cell_length_c']
        self.alpha = section['castep_cell_angle_alpha']
        self.beta  = section['castep_cell_angle_beta']
        self.gamma = section['castep_cell_angle_gamma']
        self.volume = np.sqrt( 1 - math.cos(np.deg2rad(self.alpha[0]))**2
                                 - math.cos(np.deg2rad(self.beta[0]))**2
                                 - math.cos(np.deg2rad(self.gamma[0]))**2
                                 + 2 * math.cos(np.deg2rad(self.alpha[0]))
                                     * math.cos(np.deg2rad(self.beta[0]))
                                     * math.cos(np.deg2rad(self.gamma[0])) ) * self.a[0]*self.b[0]*self.c[0]


# Storing the total energy of each SCF iteration in an array
    def onClose_section_scf_iteration(self, backend, gIndex, section):
        """trigger called when _section_scf_iteration is closed"""
        # get cached values for energy_total_scf_iteration
        ev = section['energy_total_scf_iteration']
        self.scfIterNr = len(ev)
        self.energy_total_scf_iteration_list.append(ev)

        backend.addArrayValues('energy_total_scf_iteration_list', np.asarray(self.energy_total_scf_iteration_list))
        backend.addValue('scf_dft_number_of_iterations', self.scfIterNr)


# Processing forces acting on atoms (final converged forces)
    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        #get cached values of castep_store_atom_forces
        f_st = section['castep_store_atom_forces']
        for i in range(0, self.at_nr):
            f_st[i] = f_st[i].split()
            f_st[i] = [float(j) for j in f_st[i]]
            f_st_int = f_st[i]
            self.atom_forces.append(f_st_int)
        backend.addArrayValues('atom_forces', np.asarray(self.atom_forces))
  
       
# Add SCF k points and eigenvalue from *.band file to the backend (ONLY FOR SINGLE POINT CALCULATIONS AT THIS STAGE)
        if len(self.e_spin_1) != 0:
            gIndexGroup = backend.openSection('section_eigenvalues_group')

            backend.openSection('section_eigenvalues') # opening first section_eigenvalues

            backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.k_points_scf))
            backend.addArrayValues('eigenvalues_eigenvalues', np.asarray(self.e_spin_1))
            backend.addValue('number_of_eigenvalues_kpoints', self.k_nr_scf)
            backend.addValue('number_of_eigenvalues', self.e_nr_scf)

            backend.closeSection('section_eigenvalues', gIndex)

            if len(self.e_spin_2) != 0:

                backend.openSection('section_eigenvalues') # opening the second section_eigenvalues (only for spin polarised calculations)

                backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.k_points_scf))
                backend.addArrayValues('eigenvalues_eigenvalues', np.asarray(self.e_spin_2))
                backend.addValue('number_of_eigenvalues_kpoints', self.k_nr_scf)
                backend.addValue('number_of_eigenvalues', self.e_nr_scf)

                backend.closeSection('section_eigenvalues', gIndex+1)

                backend.closeSection('section_eigenvalues_group', gIndexGroup)

            else:
                backend.closeSection('section_eigenvalues_group', gIndexGroup)
        
        backend.openSection('section_stress_tensor')
        backend.addArrayValues('stress_tensor_value',np.asarray(self.stress_tensor_value))
        backend.closeSection('section_stress_tensor', gIndex)
        
        

    
        

# Recover SCF k points and eigenvalue from *.band file (ONLY FOR SINGLE POINT CALCULATIONS AT THIS STAGE)
    def onClose_castep_section_collect_scf_eigenvalues(self, backend, gIndex, section):

        bandSuperContext = CastepBandParser.CastepBandParserContext(False)
        bandParser = AncillaryParser(
            fileDescription = CastepBandParser.build_CastepBandFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = CastepBandParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = bandSuperContext)
        
        extFile = ".bands"       # ".band_sp" = spin polarised, ".band" = not spin polarised (ONLY FOR TEST FILES IN /test/examples)
        dirName = os.path.dirname(os.path.abspath(self.fName))
        bFile = str()
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                bFile = file
                fName = os.path.normpath(os.path.join(dirName, bFile))

                with open(fName) as fIn:    
                    bandParser.parseFile(fIn)  # parsing *.band file to get SCF eigenvalues and relative k points


                self.k_nr_scf     = bandSuperContext.k_nr
                self.e_nr_scf     = bandSuperContext.e_nr
                self.k_points_scf = bandSuperContext.eigenvalues_kpoints
                self.e_spin_1     = bandSuperContext.e_spin_1
                self.e_spin_2     = bandSuperContext.e_spin_2
            else: 
                pass # if no .bands is found in the same folder as .castep file than skip and continue     
        #self.n_spin = bandSuperContext.n_spin

    def onClose_castep_section_stress_tensor(self, backend, gIndex, section): 
     #get cached values for stress tensor
        stress_tens = section['castep_store_stress_tensor']
        for i in range(len(stress_tens)):
            stress_tens[i] = stress_tens[i].split()
            stress_tens[i] = [float(j) for j in stress_tens[i]]
            stress_tens_int = stress_tens[i]
            stress_tens_int = [x / 10e9 for x in stress_tens_int] #converting GPa in Pa.
            self.stress_tensor_value.append(stress_tens_int)
    
    def onClose_castep_section_time(self, backend, gIndex, section): 
        initial_time = section['castep_initialisation_time']  
        calc_time = section['castep_calculation_time']
        final_time = section['castep_finalisation_time']
        for i in range(len(initial_time)):
            self.time_calc = (initial_time[i] + calc_time[i] + final_time[i])
        backend.addValue('time_calculation', self.time_calc)
######################################################################################
################ Triggers on closure section_system_description ######################
######################################################################################

    def onClose_section_system_description(self, backend, gIndex, section):
        """trigger called when _section_system_description is closed"""

# Processing the atom positions in fractionary coordinates (as given in the CASTEP output)
        #get cached values of castep_store_atom_position
        pos = section['castep_store_atom_position']
        self.at_nr = len(pos)
        for i in range(0, self.at_nr):
            pos[i] = pos[i].split()
            pos[i] = [float(j) for j in pos[i]]
            self.castep_atom_position.append(pos[i])
        backend.addArrayValues('castep_atom_position', np.asarray(self.castep_atom_position))


# Backend add the total number of atoms in the simulation cell
        backend.addValue('number_of_atoms', self.at_nr)


# Processing the atom labels
        #get cached values of castep_store_atom_label
        lab = section['castep_store_atom_label']
        for i in range(0, self.at_nr):
            lab[i] = re.sub('\s+', ' ', lab[i]).strip()
        self.atom_label.append(lab)
        backend.addArrayValues('atom_label', np.asarray(self.atom_label))


# Converting the fractional atomic positions (x) to cartesian coordinates (X) ( X = M^-1 x )
        for i in range(0, self.at_nr):

            pos_a = [   self.a[0] * self.castep_atom_position[i][0]
                      + self.b[0] * math.cos(np.deg2rad(self.gamma[0])) * self.castep_atom_position[i][1]
                      + self.c[0] * math.cos(np.deg2rad(self.beta[0])) * self.castep_atom_position[i][2],

                        self.b[0] * math.sin(self.gamma[0]) * self.castep_atom_position[i][1]
                      + self.c[0] * self.castep_atom_position[i][2] * (( math.cos(np.deg2rad(self.alpha[0]))
                      - math.cos(np.deg2rad(self.beta[0])) * math.cos(np.deg2rad(self.gamma[0])) ) / math.sin(np.deg2rad(self.gamma[0])) ),

                       (self.volume / (self.a[0]*self.b[0] * math.sin(np.deg2rad(self.gamma[0])))) * self.castep_atom_position[i][2] ]

            self.atom_position.append(pos_a)
        backend.addArrayValues('atom_position', np.asarray(self.atom_position), unit='angstrom')


# Backend add the simulation cell
        backend.addArrayValues('simulation_cell', np.asarray(self.cell), unit='angstrom')


######################################################################################
###################### Storing k band points and band energies #######################
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
###################### Storing k band points and band energies #######################
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
########### BAND STRUCTURE ###########################################################
######################################################################################

    def onClose_section_k_band(self, backend, gIndex, section):
        """Trigger called when section_k_band is closed.

           The k path coordinates are extracted from the *.cell input file.
        """

        cellSuperContext = CastepCellParser.CastepCellParserContext(False)
        cellParser = AncillaryParser(
            fileDescription = CastepCellParser.build_CastepCellFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = CastepCellParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = cellSuperContext)

        extFile = ".cell"       # Find the file with extension .cell
        dirName = os.path.dirname(os.path.abspath(self.fName))
        cFile = str()
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                cFile = file
        fName = os.path.normpath(os.path.join(dirName, cFile))

        with open(fName) as fIn:
            cellParser.parseFile(fIn)  # parsing *.cell file to get the k path segments
           
        self.k_start_end = cellSuperContext.k_sgt_start_end  # recover k path segments coordinartes from *.cell file
        self.k_path_nr = len(self.k_start_end)
        
        ########################################################################################
        def get_last_index(el, check):  # function that returns end index for each k path segment
            found = None
            for i, next in enumerate(check):
                if next == el:
                    found = i

            assert found != None
            return found
        ########################################################################################
        if self.k_start_end: 
            if self.castep_band_energies_1 != []:  # handling k band energies
                for i in range(self.k_nr):
                    a = [ self.castep_band_energies[i], self.castep_band_energies_1[i] ]  # spin polarised
                    self.band_en.append(a)
            else:
                self.band_en = self.castep_band_energies  # single spin

    
            path_end_index = []
            for i in range(self.k_path_nr):
                boundary = self.k_start_end[i][1]
                a = get_last_index(boundary, self.castep_band_kpoints)
                path_end_index.append(a)

            path_end_index = [0] + path_end_index  # list storing the end index of each k segment


            k_point_path = []
            for i in range(self.k_path_nr):
                a = self.castep_band_kpoints[ path_end_index[i] : path_end_index[i+1]+1 ]
                k_point_path.append(a)          # storing the k point fractional coordinates for each segment


            band_en_path = []
            for i in range(self.k_path_nr):
                a = self.band_en[ path_end_index[i] : path_end_index[i+1]+1 ]
                band_en_path.append(a)          # storing the band energies for each segment, k point and spin channel


            backend.addArrayValues('band_k_points', np.asarray(k_point_path))
            backend.addArrayValues('band_energies', np.asarray(band_en_path))
        else: 
            pass    

    
################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
######################  MAIN PARSER STARTS HERE  ###############################################################################################################
################################################################################################################################################################
############################################################################ CASTEP.Parser Version 1.0 #########################################################
################################################################################################################################################################


def build_CastepMainFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the main file of CASTEP.

    First, several subMatchers are defined, which are then used to piece together
    the final SimpleMatcher.

    Returns:
       SimpleMatcher that parses main file of CASTEP.
    """

    ########################################
  
    ########################################
    # submatcher for
    scfEigenvaluesSubMatcher = SM(name = 'scfEigenvalues',
       startReStr = r"\stype\sof\scalculation\s*\:\ssingle\spoint\senergy\s*",
       sections = ["castep_section_collect_scf_eigenvalues"]

       ) # CLOSING SM scfEigenvalues


    ########################################
    # submatcher for section method
    calculationMethodSubMatcher = SM(name = 'calculationMethods',
        startReStr = r"\susing functional\s*\:",
        forwardMatch = True,
        sections = ["section_method"],
        subMatchers = [

           SM(name = "castepXC",
              startReStr = r"\susing functional\s*\:",
              forwardMatch = True,
              sections = ["castep_section_functionals"],
              subMatchers = [

                 SM(r"\susing functional\s*\: *(?P<castep_functional_name> [A-Za-z0-9() ]*)"),
                 SM(r"\srelativistic treatment\s*\: *(?P<castep_relativity_treatment_scf> [A-Za-z0-9() -]*)")

                             ]), # CLOSING castep_section_functionals


                      ]) # CLOSING SM calculationMethods



    ########################################
    # submatcher for section_basis_set_cell_associated

    basisSetCellAssociatedSubMatcher = SM(name = 'planeWaveBasisSet',
        startReStr = r"\sbasis set accuracy\s*",
        forwardMatch = True,
        sections = ["section_basis_set_cell_associated"],
        subMatchers = [

           SM(r"\splane wave basis set cut\-off\s*\:\s*(?P<castep_basis_set_plan_wave_cutoff>[0-9.]+)")

                      ]) # CLOSING SM planeWaveBasisSet



    ########################################
    # submatcher for section system description
    systemDescriptionSubMatcher = SM(name = "systemDescription",
        startReStr = r"\s*Unit Cell\s*",
        forwardMatch = True,
        sections = ["section_system_description"],
        subMatchers = [

           # cell information
           SM(name = 'cellInformation',
              startReStr = r"\s*Unit Cell\s*",
              forwardMatch = True,
              sections = ["castep_section_cell"],
              subMatchers = [
                  SM(r"\s*(?P<castep_cell_vector>[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s+[-+0-9.eEd]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                 #SM(r"\s*(?P<castep_cell_vector>[\d\.]+\s+[\d\.]+\s+[\d\.]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                    endReStr = "\n",
                    repeats = True),

                             ]), # CLOSING castep_section_cell


           # atomic positions and cell dimesions
           SM(startReStr = r"\s*Lattice parameters",
              forwardMatch = True,
              sections = ["castep_section_atom_position"],
              subMatchers = [

                 SM(r"\s*a \=\s*(?P<castep_cell_length_a>[\d\.]+)\s*alpha \=\s*(?P<castep_cell_angle_alpha>[\d\.]+)"),
                 SM(r"\s*b \=\s*(?P<castep_cell_length_b>[\d\.]+)\s*beta  \=\s*(?P<castep_cell_angle_beta>[\d\.]+)"),
                 SM(r"\s*c \=\s*(?P<castep_cell_length_c>[\d\.]+)\s*gamma \=\s*(?P<castep_cell_angle_gamma>[\d\.]+)"),
                 SM(r"\s*x\s*(?P<castep_store_atom_label>[A-Za-z0-9]+\s+[\d\.]+)\s*[0-9]\s*(?P<castep_store_atom_position>[\d\.]+\s+[\d\.]+\s+[\d\.]+)",
                    endReStr = "\n",
                    repeats = True)

                             ]), # CLOSING castep_section_atom_position

                      ]) # CLOSING SM systemDescription



    ########################################
    # submatcher for band structure
    bandStructureSubMatcher = SM (name = 'BandStructure',
        startReStr = r"\s*\+\s*B A N D   S T R U C T U R E   C A L C U L A T I O N\s*",
        #startReStr = r"\stype\sof\scalculation\s*\:\sband\sstructure\s*",
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
                                 SM(r"\s*\+\s*[0-9]+\s*(?P<castep_store_eigenvalues>\s+[-+0-9.eEd]+)",
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
                                 SM(r"\s*\+\s*[0-9]+\s*(?P<castep_store_eigenvalues_1>\s+[-+0-9.eEd]+)",
                                    repeats = True)

                                             ]), # CLOSING castep_section_eigenvalues_1


                                    ]), # CLOSING castep_section_k_points_1


                              ]), # CLOSING 2nd section_eigenvalues


        ]) # CLOSING SM BandStructure



    ########################################
    # return main Parser ###################
    ########################################
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

    

               SM(name = 'ProgramHeader',
                  startReStr = r"\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*",
                  subMatchers = [

                     SM(r"\s\|\sWelcome to Academic Release\s(?P<program_name>[a-zA-Z]+)* version *(?P<program_version>[0-9a-zA-Z_.]*)"),
                     SM(r"\sCompiled for *(?P<program_compilation_host>[-a-zA-Z0-9._]*)\son\s(?P<castep_program_compilation_date>[a-zA-Z,\s0-9]*)\s *(?P<castep_program_compilation_time>[0-9:]*)"),
                     SM(r"\sCompiler\: *(?P<castep_compiler>[a-zA-Z\s0-9.]*)"),
                     SM(r"\sMATHLIBS\: *(?P<castep_maths_library>[a-zA-Z0-9.() ]*)\s*"),
                     SM(r"\sFFT Lib \: *(?P<castep_fft_library>[a-zA-Z0-9.() ]*)\s*"),
                     SM(r"\sFundamental constants values\: *(?P<castep_constants_reference>[a-zA-Z0-9.() ]*)\s*"),

                                  ]), # CLOSING SM ProgramHeader
                  
               scfEigenvaluesSubMatcher, # export section_eigenvalues_group to the correct nesting


               calculationMethodSubMatcher, # section_method


               basisSetCellAssociatedSubMatcher, # section_basis_set_cell_associated


               systemDescriptionSubMatcher, # section_system_description subMatcher
               SM(name = 'Atom_topology',
                  startReStr = r"\s*Mass of species in AMU\s*",              
                  endReStr = "\n",
                  #forwardMatch = True,
                  sections = ['section_topology'],
                  subMatchers = [
                 
                      SM(r"\s*(?P<castep_store_atom_name>[a-zA-Z]+)\s*(?P<castep_store_atom_mass>[\d\.]+)\s*",
                        #endReStr = "\n",   
                        repeats = True),
                     ]), # CLOSING section_atom_topology
               
                SM(startReStr = r"\-*\s*\<\-\-\sSCF\s*",
                  forwardMatch = True,
                  sections = ["section_single_configuration_calculation"],
                  subMatchers = [                     
                 
                     SM(name = 'ScfIterations',
                        startReStr = r"SCF\sloop\s*Energy\s*Fermi\s*Energy\sgain\s*Timer\s*<\-\-\sSCF\s*",
                        sections = ['section_scf_iteration'],
                        subMatchers = [
                     
                           SM(r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                              repeats = True),

                                  ]), # CLOSING section_scf_iteration

                     SM(name = 'ScfIterations',
                        startReStr = r"SCF\sloop\s*Energy\s*Energy\sgain\s*Timer\s*<\-\-\sSCF\s*",
                        sections = ['section_scf_iteration'],
                        subMatchers = [

                           SM(r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                              repeats = True),

                                  ]), # CLOSING section_scf_iteration

                      SM(r"Final energy = *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                      SM(r"Final energy\,\s*E\s*= *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                      SM(r"Final free energy\s*\(E\-TS\)\s*= *(?P<energy_free__eV>[-+0-9.eEdD]*)"), # matching final converged total free energy
                      SM(r"NB est\. 0K energy\s*\(E\-0\.5TS\)\s*= *(?P<energy_total_T0__eV>[-+0-9.eEdD]*)"), # 0K corrected final SCF energy

                     bandStructureSubMatcher,  # band structure subMatcher
                     
                     

                     SM(startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* Forces \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\s*",
                        subMatchers = [
                           SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True)
                                      ]),

                     SM(startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* Symmetrised Forces \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\s*",
                        subMatchers = [
                           SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True)
                                      ]),
                     
                     SM(name = 'stresstensor',
                        startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\* Symmetrised Stress Tensor \*\*\*\*\*\*\*\*\*\*\*\s*",
                        sections = ['castep_section_stress_tensor'],
                        subMatchers = [
                           SM(r"\s*\*\s*[a-z]\s*(?P<castep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True),  
                                      ]), # CLOSING section_stress_tensor
                     
                     SM(name = 'calc_time',
                        startReStr = r" A BibTeX formatted list of references used in this run has been written to",
                        sections = ['castep_section_time'],
                        subMatchers = [
                           SM(r"Initialisation time =\s*(?P<castep_initialisation_time>[0-9.]*)"), # matching final converged total energy
                           SM(r"Calculation time    =\s*(?P<castep_calculation_time>[0-9.]*)"),
                           SM(r"Finalisation time   =\s*(?P<castep_finalisation_time>[0-9.]*)"),
                                      ]), # CLOSING section_castep_time
                                        ]) # CLOSING section_single_configuration_calculation                
                
                
                           ]) # CLOSING SM NewRun        
                
 
        ])



def get_cachingLevelForMetaName(metaInfoEnv):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
                                #'band_energies' : CachingLevel.Cache,
                                #'band_k_points' : CachingLevel.Cache,
                                'castep_basis_set_plan_wave_cutoff' : CachingLevel.Cache,
                                #'eigenvalues_eigenvalues': CachingLevel.Cache,
                                'eigenvalues_kpoints':CachingLevel.Cache
                                }

    # Set caching for temparary storage variables
    for name in metaInfoEnv.infoKinds:
        if (   name.startswith('castep_store_')
            or name.startswith('castep_cell_')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main():
    """Main function.

    Set up everything for the parsing of the CASTEP main file and run the parsing.
    """
    # get main file description
    CastepMainFileSimpleMatcher = build_CastepMainFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json
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
                 defaultSectionCachingLevel = True)

if __name__ == "__main__":
    main()








