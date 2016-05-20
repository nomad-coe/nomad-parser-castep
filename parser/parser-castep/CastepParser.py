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
import CastepMDParser
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
        self.func_total                        = []
        self.relativistic                      = []
        self.dispersion                        = None
        self.cell                              = []
        self.at_nr                             = 0
        self.atom_type_mass                    = []
        self.atom_labels                       = []
        self.atom_optim_position = []
        self.castep_optimised_atom_positions   = [] 
        self.castep_atom_optim_label           = []
        self.atom_forces                       = []
        self.atom_forces_band                       = []
        self.stress_tensor_value               = []
        self.raman_tensor_value               = []
        self.castep_atom_positions              = []
        self.atom_positions                     = []
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
        self.functional_type                   = []
        self.functional_weight                 = None
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
        self.van_der_waals_name                = None
        self.e_spin_1                          = []
        self.e_spin_2                          = []
        self.optim_method                      = []
        self.disp_energy                       = [] 
        self.geoConvergence = None
        self.total_energy_corrected_for_finite_basis = []
        self.contr_s =[]
        self.contr_p =[]
        self.contr_d = []
        self.contr_f = []
        self.total_contribution = []
        self.total_charge = []
        self.frequencies = []
        self.ir_intens = []
        self.raman_act = []
        self.irr_repres = []
        self.nr_iter =[]
        self.castep_atom_velocities = []
        self.time_0 = []
        self.frame_temp = []
        self.frame_press =[]
        self.frame_total =[]
        self.frame_total_energy = []
        self.frame_hamiltonian = []
        self.frame_kinetic = []
        self.frame_atom_forces=[]
        self.frame_atom_lables =[]
        self.frame_stress_tensor =[]
        self.frame_position =[]
        self.frame_cell=[]
        self.frame_time=[]
        self.energy_frame =[]
        self.total_energy_frame = []
        self.frame_energies =[]
        self.scfgIndex=[]
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
        
        #relativistic_names = section["castep_relativity_treatment_scf"]

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
            " LDA-X" : "LDA_X_PZ",
            " LDA-C" : "LDA_C_PZ",
            #" 'EXX','EXX-LDA'" : ""
            #" 'SHF','SX'"      :
            #" Optimised Effective Potential" :
        }
       
       
        # Match each castep functional name and sort the matches into a list
        self.functionals = []

        for name in functional_names:
            match = functional_map.get(name)
            if match:
                self.functionals.append(match)
        #self.functionals = "_".join(sorted(self.functionals))
        

    def onClose_castep_section_functional_definition(self, backend, gIndex, section):

        self.functional_types = section["castep_functional_type"]
        self.functional_weight = section["castep_functional_weight"]
        

        # Define a mapping for the functionals
        functional_map = {
            "PBE": "GGA_C_PBE_GGA_X_PBE",
            "PW91": "GGA_C_PW91_GGA_X_PW91",
            "RPBE": "GGA_X_RPBE",
            "PBEsol": "GGA_X_PBE_SOL",          
            "HF"       :"HF_X", 
            "HF-LDA":"HF_X_LDA_C_PW",
            "HSE03"              :"HYB_GGA_XC_HSE03",
            "HSE06"              :"HYB_GGA_XC_HSE06",
            "PBE0"               :"GGA_C_PBE_GGA_X_PBE", 
            "WC" :"GGA_C_PBE_GGA_X_WC",    
            "LDA-X" : "LDA_X_PZ",
            "LDA-C" : "LDA_C_PZ",
            "LDA": "LDA_C_PZ_LDA_X_PZ",
            "b3lyp"       :"HYB_GGA_XC_B3LYP5",
            #" 'EXX','EXX-LDA'" : ""
            #" 'SHF','SX'"      :
            #" Optimised Effective Potential" :
        }
        
        self.functionals = []
       
        for name in self.functional_types:
            match = functional_map.get(name)
            if match:
                self.functionals.append(match)
        #self.functionals = "_".join(sorted(self.functionals))
        
    def onClose_castep_section_relativity_treatment(self, backend, gIndex, section):
        relativistic_names = section["castep_relativity_treatment_scf"]
        
        # Define a mapping for the relativistic treatments
        relativistic_map = {
            " Koelling-Harmon": "scalar_relativistic"
        }
        self.relativistic = []

        for name in relativistic_names:
            match = relativistic_map.get(name)
            if match:
                self.relativistic.append(match)
        self.relativistic = "_".join(sorted(self.relativistic))
# Here we add info about the XC functional and relativistic treatment
    
    def onClose_castep_section_geom_optimisation_method(self, backend, gIndex, section):
        self.optim_method = section["castep_geometry_optim_method"]

    def onClose_section_method(self, backend, gIndex, section):
        self.van_der_waals_name = section["van_der_Waals_method"]
        
        if self.van_der_waals_name is not None:
            dispersion_map = {
            " OBS correction scheme": "OBS",
            " G06 correction scheme": "G06",
            " TS correction scheme" : "TS",
            " JCHS correction scheme" : "JCHS",
            }
            self.dispersion = []
            for name in self.van_der_waals_name:
                match = dispersion_map.get(name)
                if match:
                    self.dispersion.append(match)    
            self.dispersion = "_".join(sorted(self.dispersion))
        else:
            pass
        #backend.addValue('van_der_Waals_method',self.dispersion)
        
        
        if self.functional_weight is not None:
            self.func_and_weight = []
            for i in range(len(self.functional_types)):
                self.func_total.append(self.functionals[i]+'_'+self.functional_weight[i]) 
                backend.openSection('section_XC_functionals')            
                backend.addValue('XC_functional_name', self.functionals[i]) 
                backend.addValue('XC_functional_weight', self.functional_weight[i])
                backend.closeSection('section_XC_functionals',gIndex+i)        
        # Push the functional string into the backend
        # Push the relativistic treatment string into the backend
            backend.addValue('XC_functional', "_".join(sorted(self.functionals)))
            backend.addValue('relativity_method', self.relativistic)
        #if self.functional_weight = 0
            if self.dispersion is not None:
                backend.addValue('XC_method_current', ("_".join(sorted(self.func_total)))+'_'+self.dispersion+'_'+self.relativistic)
            else:
                backend.addValue('XC_method_current', ("_".join(sorted(self.func_total)))+'_'+self.relativistic)
        else:
            for i in range(len(self.functionals)):
          #      self.func_total.append(self.functionals[i]+'_'+self.functional_weight[i])
                backend.openSection('section_XC_functionals')            
                backend.addValue('XC_functional_name', self.functionals[i]) 
         #       backend.addValue('XC_functional_weight', self.functional_weight[i])
                backend.closeSection('section_XC_functionals',gIndex+i)
            backend.addValue('XC_functional', "_".join(sorted(self.functionals)))
            backend.addValue('relativity_method', self.relativistic)
            if self.dispersion is not None:
                backend.addValue('XC_method_current', ("_".join(sorted(self.functionals)))+'_'+self.dispersion+'_'+self.relativistic)
            else:
                backend.addValue('XC_method_current', ("_".join(sorted(self.functionals)))+'_'+self.relativistic)
               
        
# Here we add basis set name and kind for the plane wave code
    def onClose_section_basis_set_cell_dependent(self, backend, gIndex, section):
        ecut_str = section['castep_basis_set_planewave_cutoff']
        self.ecut = float(ecut_str[0])
        eVtoRy = 0.073498618
        ecut_str_name = int(round(eVtoRy*self.ecut))

        basis_set_kind = 'plane_waves'
        basis_set_name = 'PW_'+str(ecut_str_name)
        backend.addValue('basis_set_planewave_cutoff', self.ecut)
        backend.addValue('basis_set_cell_dependent_kind', basis_set_kind)
        backend.addValue('basis_set_cell_dependent_name', basis_set_name)

    def onClose_castep_section_cell_optim(self, backend, gIndex, section):
        """trigger called when _castep_section_cell is closed"""
        # get cached values for castep_cell_vector
        vet = section['castep_cell_vector_optim']
     
        vet[0] = vet[0].split()
        vet[0] = [float(i) for i in vet[0]]

        vet[1] = vet[1].split()
        vet[1] = [float(i) for i in vet[1]]

        vet[2] = vet[2].split()
        vet[2] = [float(i) for i in vet[2]]

        self.cell.append(vet[0])
        self.cell.append(vet[1])
        self.cell.append(vet[2]) # Reconstructing the unit cell vector by vector    
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
    def onClose_castep_section_atom_positions(self, backend, gIndex, section):
        """trigger called when _castep_section_atom_positions is closed"""
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

    def onClose_castep_section_atom_positions_optim(self, backend, gIndex, section):
        """trigger called when _castep_section_atom_positions is closed"""
        # get cached values for cell magnitudes and angles
    
        self.a = section['castep_cell_length_a_optim']
        self.b = section['castep_cell_length_b_optim']
        self.c = section['castep_cell_length_c_optim']
        self.alpha = section['castep_cell_angle_alpha_optim']
        self.beta  = section['castep_cell_angle_beta_optim']
        self.gamma = section['castep_cell_angle_gamma_optim']
        self.volume = np.sqrt( 1 - math.cos(np.deg2rad(self.alpha[0]))**2
                                 - math.cos(np.deg2rad(self.beta[0]))**2
                                 - math.cos(np.deg2rad(self.gamma[0]))**2
                                 + 2 * math.cos(np.deg2rad(self.alpha[0]))
                                     * math.cos(np.deg2rad(self.beta[0]))
                                     * math.cos(np.deg2rad(self.gamma[0])) ) * self.a[0]*self.b[0]*self.c[0]    
       
# Storing the total energy of each SCF iteration in an array
    # def onClose_section_scf_iteration(self, backend, gIndex, section):
    #     """trigger called when _section_scf_iteration is closed"""
    #     # get cached values for energy_total_scf_iteration
    #     ev = section['energy_total_scf_iteration']
    #     self.scfIterNr = len(ev)
    #     self.energy_total_scf_iteration_list.append(ev)

    #     backend.addArrayValues('energy_total_scf_iteration_list', np.asarray(self.energy_total_scf_iteration_list))
    #     backend.addValue('number_of_scf_iterations', self.scfIterNr)

 
# Processing forces acting on atoms (final converged forces)
    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        self.time_0 = section['castep_frame_time_0'] 

        f_st = section['castep_store_atom_forces']
      
        if f_st is not None:
            for i in range(0, self.at_nr):
                f_st[i] = f_st[i].split()
                f_st[i] = [float(j) for j in f_st[i]]
                f_st_int = f_st[i]
                 
                self.atom_forces.append(f_st_int)
                self.atom_forces = self.atom_forces[-self.at_nr:] 
                
            backend.addArrayValues('atom_forces', np.asarray(self.atom_forces))
        else: 
            pass        
# Add SCF k points and eigenvalue from *.band file to the backend (ONLY FOR SINGLE POINT CALCULATIONS AT THIS STAGE)
        if len(self.e_spin_1) != 0:
            # gIndexGroup = backend.openSection('section_eigenvalues_group')

            backend.openSection('section_eigenvalues') # opening first section_eigenvalues

            backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.k_points_scf))
            backend.addArrayValues('eigenvalues_values', np.asarray(self.e_spin_1))
            backend.addValue('number_of_eigenvalues_kpoints', self.k_nr_scf)
            backend.addValue('number_of_eigenvalues', self.e_nr_scf)

            backend.closeSection('section_eigenvalues', gIndex)

        if len(self.e_spin_2) != 0:

            backend.openSection('section_eigenvalues') # opening the second section_eigenvalues (only for spin polarised calculations)

            backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.k_points_scf))
            backend.addArrayValues('eigenvalues_values', np.asarray(self.e_spin_2))
            backend.addValue('number_of_eigenvalues_kpoints', self.k_nr_scf)
            backend.addValue('number_of_eigenvalues', self.e_nr_scf)

            backend.closeSection('section_eigenvalues', gIndex+1)

            #     backend.closeSection('section_eigenvalues_group', gIndexGroup)

            # else:
            #     backend.closeSection('section_eigenvalues_group', gIndexGroup)
        
        
        
        
        #backend.openSection('section_stress_tensor')
        if self.stress_tensor_value:
            backend.addArrayValues('stress_tensor',np.asarray(self.stress_tensor_value))
        #backend.closeSection('section_stress_tensor', gIndex)
        #backend.addValue('time_calculation', self.time_calc)

        #######Total Dispersion energy recovered from totatl energy. ##############      
        if self.dispersion is not None:
            van_der_waals_energy = section['castep_total_dispersion_corrected_energy']
            total_energy = section['energy_total']
            for i in range(len(total_energy)):
                self.disp_energy = abs(van_der_waals_energy[i] - total_energy[i])
            backend.addValue('energy_van_der_Waals', self.disp_energy)

        
    def onClose_castep_section_SCF_iteration_frame(self, backend, gIndex, section):
        self.frame_energies = section['castep_SCF_frame_energy']
                       
        frame_time = section['castep_frame_time']

        if self.frame_energies:
                
                # self.energy_frame =[]
            for i in range(len(self.frame_energies)):
                     
                self.frame_energies[i]=self.frame_energies[i].split()
                self.frame_energies[i]=[float(j) for j in self.frame_energies[i]]
                energies = self.frame_energies[i]
                self.energy_frame.append(energies)   

            if frame_time:
                for i in range(len(frame_time)):                   
                    self.time_0.extend(frame_time)                
                    
        #get cached values of castep_store_atom_forces
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
        stress_tens =[]
        stress_tens = section['castep_store_stress_tensor']
    
        for i in range(len(stress_tens)):
            stress_tens[i] = stress_tens[i].split()
            stress_tens[i] = [float(j) for j in stress_tens[i]]
            stress_tens_int = stress_tens[i]
            stress_tens_int = [x / 10e9 for x in stress_tens_int] #converting GPa in Pa.
            self.stress_tensor_value.append(stress_tens_int)
        self.stress_tensor_value = self.stress_tensor_value[-3:]     

    def onClose_section_castep_raman_tensor(self, backend, gIndex, section): 
     #get cached values for stress tensor
        ram_tens =[]
        ram_tens = section['castep_store_raman_tensor']
        
        for i in range(len(ram_tens)):
            ram_tens[i] = ram_tens[i].split()
            ram_tens[i] = [float(j) for j in ram_tens[i]]
            ram_tens_int = ram_tens[i]
            self.raman_tensor_value.append(ram_tens_int)
        self.raman_tensor_value = self.raman_tensor_value[-3:]         
        if self.raman_tensor_value:
            backend.addArrayValues('castep_ramen_tensor',np.asarray(self.raman_tensor_value))    
    
    def onClose_castep_section_time(self, backend, gIndex, section): 
        initial_time = section['castep_initialisation_time']  
        calc_time = section['castep_calculation_time']
        final_time = section['castep_finalisation_time']
        for i in range(len(initial_time)):
            self.time_calc = (initial_time[i] + calc_time[i] + final_time[i])
            

######################################################################################
################ Triggers on closure section_system ######################
######################################################################################

    def onClose_section_system(self, backend, gIndex, section):
        """trigger called when _section_system is closed"""

# Processing the atom positions in fractionary coordinates (as given in the CASTEP output)
        #get cached values of castep_store_atom_positions
        n_electrons = section['castep_number_of_electrons_store']
        if n_electrons:
            backend.addValue('castep_number_of_electrons', n_electrons)
        pos = section['castep_store_atom_positions']
        if pos:
            self.at_nr = len(pos)
            for i in range(0, self.at_nr):
                pos[i] = pos[i].split()
                pos[i] = [float(j) for j in pos[i]]
                self.castep_atom_positions.append(pos[i])
            backend.addArrayValues('castep_atom_positions', np.asarray(self.castep_atom_positions))


# Backend add the total number of atoms in the simulation cell
            backend.addValue('number_of_atoms', self.at_nr)


# Processing the atom labels
        #get cached values of castep_store_atom_labels
            lab = section['castep_store_atom_labels']
            for i in range(0, self.at_nr):
                lab[i] = re.sub('\s+', ' ', lab[i]).strip()
            self.atom_labels.append(lab)
            backend.addArrayValues('atom_labels', np.asarray(self.atom_labels))


# Converting the fractional atomic positions (x) to cartesian coordinates (X) ( X = M^-1 x )
            for i in range(0, self.at_nr):

                pos_a = [   self.a[0] * self.castep_atom_positions[i][0]
                      + self.b[0] * math.cos(np.deg2rad(self.gamma[0])) * self.castep_atom_positions[i][1]
                      + self.c[0] * math.cos(np.deg2rad(self.beta[0])) * self.castep_atom_positions[i][2],

                        self.b[0] * math.sin(self.gamma[0]) * self.castep_atom_positions[i][1]
                      + self.c[0] * self.castep_atom_positions[i][2] * (( math.cos(np.deg2rad(self.alpha[0]))
                      - math.cos(np.deg2rad(self.beta[0])) * math.cos(np.deg2rad(self.gamma[0])) ) / math.sin(np.deg2rad(self.gamma[0])) ),

                       (self.volume / (self.a[0]*self.b[0] * math.sin(np.deg2rad(self.gamma[0])))) * self.castep_atom_positions[i][2] ]

                self.atom_positions.append(pos_a)
            
            backend.addArrayValues('simulation_cell', np.asarray(self.cell[-3:]))
            backend.addValue('CASTEP_cell_volume', self.volume)     
            backend.addArrayValues('atom_positions', np.asarray(self.atom_positions))
        
        vel = section['castep_store_atom_ionic_velocities']
        if vel:
            for i in range(0, self.at_nr):
                # vel[i] = [float(j) for j in vel[i]]
                self.castep_atom_velocities.append(vel[i])
            backend.addArrayValues('castep_atom_ionic_velocities', np.asarray(self.castep_atom_velocities))

# Backend add the simulation cell
        pos_opt = section['castep_store_optimised_atom_positions']
        if pos_opt:

            backend.addArrayValues('atom_labels', np.asarray(self.atom_labels[-self.at_nr:]))
            
            self.at_nr_opt = len(pos_opt)
            for i in range(0, self.at_nr_opt):
                pos_opt[i] = pos_opt[i].split()
                pos_opt[i] = [float(j) for j in pos_opt[i]]
                self.castep_optimised_atom_positions.append(pos_opt[i])
            backend.addArrayValues('castep_atom_positions', np.asarray(self.castep_optimised_atom_positions[-self.at_nr_opt:]))
        # #     print pos_opt[i]    
       
# Converting the fractional atomic positions (x) to cartesian coordinates (X) ( X = M^-1 x )
            for i in range(0, self.at_nr_opt):

                pos_opt_a = [   self.a[0] * self.castep_optimised_atom_positions[i][0]
                      + self.b[0] * math.cos(np.deg2rad(self.gamma[0])) * self.castep_optimised_atom_positions[i][1]
                      + self.c[0] * math.cos(np.deg2rad(self.beta[0])) * self.castep_optimised_atom_positions[i][2],

                        self.b[0] * math.sin(self.gamma[0]) * self.castep_optimised_atom_positions[i][1]
                      + self.c[0] * self.castep_optimised_atom_positions[i][2] * (( math.cos(np.deg2rad(self.alpha[0]))
                      - math.cos(np.deg2rad(self.beta[0])) * math.cos(np.deg2rad(self.gamma[0])) ) / math.sin(np.deg2rad(self.gamma[0])) ),

                       (self.volume / (self.a[0]*self.b[0] * math.sin(np.deg2rad(self.gamma[0])))) * self.castep_optimised_atom_positions[i][2] ]

                self.atom_optim_position.append(pos_opt_a)
            
            backend.addArrayValues('simulation_cell', np.asarray(self.cell[-3:])) 
            backend.addArrayValues('atom_positions', np.asarray(self.atom_optim_position[-self.at_nr_opt:]))
            backend.addValue('CASTEP_cell_volume', self.volume) 
        else:
            pass
        
        # if self.volume:
        #     backend.addValue('CASTEP_cell_volume', self.volume) 
       
        # if self.cell: 
        #     backend.addArrayValues('simulation_cell', np.asarray(self.cell[-3:]))

    def onClose_castep_section_vibrational_frequencies(self, backend, gIndex, section):
        frequ = section['castep_vibrationl_frequencies_store']
        ##### in this case the name castep_ir_store sands for irreducible representation in the point group.
        irr_rep = section ['castep_ir_store']
        ##### in this case the name castep_ir_intesities_store sands for Infra-red intensities. 
        ir_intensities = section['castep_ir_intensity_store']
     
        self.nr_iter =section['castep_n_iterations_phonons']

        raman_activity = section['castep_raman_activity_store']
        
        #raman_act = section['castep_raman_active_store']
        
        if frequ and ir_intensities and raman_activity:
                for i in range(0,len(frequ)):
                    frequ[i] = frequ[i].split()
                    frequ[i] = [float(j) for j in frequ[i]]
                    frequ_list = frequ[i]
                    self.frequencies.append(frequ_list)    
                    irr_rep[i] = irr_rep[i].split()
                    irr_rep_list = irr_rep[i]
                    self.irr_repres.append(irr_rep_list)
                    ir_intensities[i] = ir_intensities[i].split()
                    ir_intensities[i] = [float(j) for j in ir_intensities[i]]
                    ir_intens_list = ir_intensities[i]
                    self.ir_intens.append(ir_intens_list) 
                    
                    raman_activity[i] = raman_activity[i].split()
                   
                    raman_activity[i] = [float(j) for j in raman_activity[i]]
                    raman_list = raman_activity[i]
                    self.raman_act.append(raman_list)  
                    
                backend.addArrayValues('castep_ir_intensity', np.asarray(self.ir_intens[-len(self.nr_iter):]))
                backend.addArrayValues('castep_vibrationl_frequencies', np.asarray(self.frequencies[-len(self.nr_iter):])) 
                backend.addArrayValues('castep_ir', np.asarray(self.irr_repres[-len(self.nr_iter):]))
                backend.addArrayValues('castep_raman_activity', np.asarray(self.raman_act[-len(self.nr_iter):]))
        elif frequ and ir_intensities:
                for i in range(0,len(frequ)):
                    frequ[i] = frequ[i].split()
                    frequ[i] = [float(j) for j in frequ[i]]
                    frequ_list = frequ[i]
                    self.frequencies.append(frequ_list)    
                    irr_rep[i] = irr_rep[i].split()
                    irr_rep_list = irr_rep[i]
                    self.irr_repres.append(irr_rep_list)
                    ir_intensities[i] = ir_intensities[i].split()
                    ir_intensities[i] = [float(j) for j in ir_intensities[i]]
                    ir_intens_list = ir_intensities[i]
                    self.ir_intens.append(ir_intens_list)
                backend.addArrayValues('castep_ir_intensity', np.asarray(self.ir_intens[-len(self.nr_iter):]))
                backend.addArrayValues('castep_vibrationl_frequencies', np.asarray(self.frequencies[-len(self.nr_iter):]))    
                backend.addArrayValues('castep_ir', np.asarray(self.irr_repres[-len(self.nr_iter):]))
        elif frequ:
                for i in range(0,len(frequ)):
                    frequ[i] = frequ[i].split()
                    frequ[i] = [float(j) for j in frequ[i]]
                    frequ_list = frequ[i]
                    self.frequencies.append(frequ_list)    
                    irr_rep[i] = irr_rep[i].split()
                    irr_rep_list = irr_rep[i]
                    self.irr_repres.append(irr_rep_list) 
                backend.addArrayValues('castep_vibrationl_frequencies', np.asarray(self.frequencies[-len(self.nr_iter):])) 
                backend.addArrayValues('castep_ir', np.asarray(self.irr_repres[-len(self.nr_iter):]))
            
    def onClose_castep_section_population_analysis(self, backend, gIndex, section):
        orb_contr = section['castep_orbital_contributions']
        tot_charge = section ['castep_mulliken_charge_store']
        
        if orb_contr:
            for i in range(0, self.at_nr):
                
                orb_contr[i] = orb_contr[i].split()
                orb_contr[i] = [float(j) for j in orb_contr[i]]
                orb_contr_a = orb_contr[i]
                tot= orb_contr_a[0]+orb_contr_a[1]+orb_contr_a[2]+orb_contr_a[3]
                
                self.total_contribution.append(tot)
                
                self.contr_s.append(orb_contr_a[0])
                self.contr_p.append(orb_contr_a[1])
                self.contr_d.append(orb_contr_a[2])
                self.contr_f.append(orb_contr_a[3])
                 
            
            backend.addArrayValues('castep_orbital_s', np.asarray(self.contr_s))
            backend.addArrayValues('castep_orbital_p', np.asarray(self.contr_p))
            backend.addArrayValues('castep_orbital_d', np.asarray(self.contr_d))
            backend.addArrayValues('castep_orbital_f', np.asarray(self.contr_f))
            backend.addValue('castep_total_orbital',self.total_contribution)
        
        if tot_charge:
            for i in range(0, self.at_nr):                
                tot_charge[i] = tot_charge[i].split()
                tot_charge[i] = [float(j) for j in tot_charge[i]]
                total_charge_list = tot_charge[i]
                self.total_charge.append(total_charge_list)
            backend.addArrayValues('castep_mulliken_charge', np.asarray(self.total_charge))

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

            for i in range(self.k_path_nr):    
                backend.openSection('section_k_band_segment')
                backend.addArrayValues('band_k_points', np.asarray(k_point_path[i]))
                backend.addArrayValues('band_energies', np.asarray(band_en_path[i]))
                backend.addArrayValues('band_segm_start_end', np.asarray(self.k_start_end[i])) 
                backend.addValue('number_of_k_points_per_segment', len(k_point_path[i])) 
                backend.closeSection('section_k_band_segment',i)
        else: 
            pass    

    def onClose_section_run(self, backend, gIndex, section):
        f_st = section['castep_store_atom_forces_band']
        time_list = self.time_0
        if f_st:
            for i in range(0, self.at_nr):
                f_st[i] = f_st[i].split()
                f_st[i] = [float(j) for j in f_st[i]]

                f_st_int = f_st[i]
                 
                self.atom_forces_band.append(f_st_int)
                self.atom_forces_band = self.atom_forces_band[-self.at_nr:] 
            backend.addArrayValues('castep_atom_forces', np.asarray(self.atom_forces_band))
        else: 
            pass        
        #optim_succ = section['CASTEP_geom_converged']
        if section['CASTEP_geom_converged'] is not None:
            if section['CASTEP_geom_converged'][-1] == 'successfully':
                self.geoConvergence = True
            else:
                self.geoConvergence = False
        if self.geoConvergence is not None:
            backend.openSection('section_frame_sequence')
            backend.addValue('geometry_optimization_converged', self.geoConvergence)        
            backend.closeSection('section_frame_sequence',gIndex)
        MDSuperContext = CastepMDParser.CastepMDParserContext(False)
        MDParser = AncillaryParser(
            fileDescription = CastepMDParser.build_CastepMDFileSimpleMatcher(),
            parser = self.parser,
            cachingLevelForMetaName = CastepMDParser.get_cachingLevelForMetaName(self.metaInfoEnv, CachingLevel.Ignore),
            superContext = MDSuperContext)

        extFile = ".md"       # Find the file with extension .cell
        dirName = os.path.dirname(os.path.abspath(self.fName))
        cFile = str()
        for file in os.listdir(dirName):
            if file.endswith(extFile):
                cFile = file
        fName = os.path.normpath(os.path.join(dirName, cFile))

         # parsing *.cell file to get the k path segments
        if file.endswith(extFile):   
            with open(fName) as fIn:
                MDParser.parseFile(fIn)
            
            self.frame_temp = MDSuperContext.frame_temperature  # recover k path segments coordinartes from *.cell file
            self.frame_press= MDSuperContext.frame_pressure
            self.frame_kinetic = MDSuperContext.kinetic
            self.frame_total = MDSuperContext.total_energy
            self.frame_atom_forces = MDSuperContext.total_forces
            self.frame_atom_veloc = MDSuperContext.total_velocities
            self.frame_stress_tensor = MDSuperContext.frame_stress_tensor
            self.frame_position =MDSuperContext.total_positions
            self.frame_cell=MDSuperContext.frame_cell
            self.frame_vet_velocities =MDSuperContext.vector_velocities
            

            
            if self.frame_temp:
                
                gIndexGroupscf = backend.openSection('section_scf_iteration')
                for i in range(len(self.frame_atom_forces)):
                   
                    backend.openSection('section_system')
                    backend.addArrayValues('atom_velocities', np.asarray(self.frame_atom_veloc[i]))
                    backend.addArrayValues('atom_labels', np.asarray(self.atom_labels))
                    backend.addArrayValues('atom_positions', np.asarray(self.frame_position[i]))
                    backend.addArrayValues('simulation_cell', np.asarray(self.frame_cell[i]))
                    backend.addArrayValues('castep_velocities_cell_vector',np.asarray(self.frame_vet_velocities[i]))
                    backend.closeSection('section_system',i+2)
                    
                    backend.openSection('section_single_configuration_calculation')
                    backend.addArrayValues('atom_forces', np.asarray(self.frame_atom_forces[i]))
                    backend.addArrayValues('stress_tensor',np.asarray(self.frame_stress_tensor[i]))
                    
                    
                    
                    if i > 0:    
                        for j in range(len(self.frame_energies)):
                            s = j + i*len(self.frame_energies) - len(self.frame_energies)
                            
                            backend.openSection('section_scf_iteration')
                            backend.addValue('energy_total_scf_iteration', self.energy_frame[s])
                            
                            backend.closeSection('section_scf_iteration',s+gIndexGroupscf)
                                     
                    
                    backend.closeSection('section_single_configuration_calculation',i+1) 

                backend.openSection('section_frame_sequence')
                backend.addValue('number_of_frames_in_sequence',(len(self.frame_total)))
                backend.addArrayValues('frame_sequence_temperature', np.asarray(self.frame_temp))
                backend.addArrayValues('frame_sequence_pressure', np.asarray(self.frame_press))
                backend.addArrayValues('frame_sequence_kinetic_energy', np.asarray(self.frame_kinetic))
                backend.addArrayValues('frame_sequence_potential_energy', np.asarray(self.frame_total))
                backend.addArrayValues('frame_sequence_time', np.asarray(time_list))
                backend.closeSection('section_frame_sequence',gIndex)
            else:
                pass
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
        startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* Exchange-Correlation Parameters \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*",
        forwardMatch = True,
        sections = ["section_method"],
        subMatchers = [

           SM(name = "castepXC",
              startReStr = r"\susing functional\s*\:",
              forwardMatch = True,
              sections = ["castep_section_functionals"],
              subMatchers = [

                 SM(r"\susing functional\s*\: *(?P<castep_functional_name> [A-Za-z0-9() ]*)"),
                 #SM(r"\srelativistic treatment\s*\: *(?P<castep_relativity_treatment_scf> [A-Za-z0-9() -]*)")
                             ]), # CLOSING castep_section_functionals
            
           SM(name = "castepXC_definition",
              startReStr = r"\susing custom XC functional definition\:",
              #endReStr = r"\srelativistic treatment\s*\:\s*",
              forwardMatch = True,
              sections = ["castep_section_functional_definition"],
              subMatchers = [     
                 SM(r"\s*(?P<castep_functional_type>[A-Za-z0-9]+)\s*(?P<castep_functional_weight>[0-9.]+)",
                        repeats = True),
                 #SM(r"\srelativistic treatment\s*\: *(?P<castep_relativity_treatment_scf> [A-Za-z0-9() -]*)")              
                              ]), # CLOSING castep_section_functional_definition

           SM(name = "castep_relativ",
              startReStr = r"\srelativistic treatment\s*\:\s*",
              forwardMatch = True,
              sections = ["castep_section_relativity_treatment"],
              subMatchers = [
                 SM(r"\srelativistic treatment\s*\: *(?P<castep_relativity_treatment_scf> [A-Za-z0-9() -]+)")
                             ]), # CLOSING castep_section_relativistic_treatment           

            SM(name = "van der Waals",
               startReStr = r"\sDFT\+D: Semi-empirical dispersion correction\s*\:\s*",
               #forwardMatch = True,
               #sections = ["castep_section_relativity_treatment"],
               subMatchers = [
                 SM(r"\sSEDC with\s*\: *(?P<van_der_Waals_method> [A-Za-z0-9() -]+)"),
                             ]), # CLOSING van der waals
            
            
              

              ]) # CLOSING SM calculationMethods



    ########################################
    # submatcher for section_basis_set_cell_dependent

    basisSetCellAssociatedSubMatcher = SM(name = 'planeWaveBasisSet',
        startReStr = r"\s\*\*\** Basis Set Parameters \*\*\**\s*",
        forwardMatch = True,
        sections = ["section_basis_set_cell_dependent"],
        subMatchers = [

            SM(r"\splane wave basis set cut\-off\s*\:\s*(?P<castep_basis_set_planewave_cutoff>[0-9.]+)")
            
                      ]) # CLOSING SM planeWaveBasisSet


    phononCalculationSubMatcher = SM(name = 'phonon_calculation',
        sections = ["castep_section_phonons"],
        startReStr = r"\s\*\*\** Phonon Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\sphonon calculation method\s*\:\s*(?P<castep_phonon_method>[a-zA-Z]+\s[a-zA-Z]+)"),
            SM(r"\sDFPT solver method\s*\:\s*(?P<castep_DFPT_solver_method>[a-zA-Z0-9.() ]+)"),
                          ])

    GeomOptimParameterSubMatcher = SM(name = 'optimistation_parameters',
        sections = ["castep_section_geom_optimisation_method"],
        startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\ Geometry Optimization Parameters \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\s*",
        subMatchers = [
            
            SM(r"\soptimization method\s*\:\s*(?P<castep_geometry_optim_method>[a-zA-Z]+)"),
                                            ]) # CLOSING converged optmisation
    
    ElectronicParameterSubMatcher = SM(name = 'Elec_parameters' ,            
        sections = ["section_system"],
        startReStr = r"\s\*\*\** Electronic Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*number of  electrons\s*\:\s*(?P<castep_number_of_electrons_store>[0-9.]+)"),
            ])
    
    MDParameterSubMatcher = SM(name = 'MD_parameters' ,            
        sections = ["section_sampling_method"],
        startReStr = r"\s\*\*\** Molecular Dynamics Parameters \*\*\**\s*",
        subMatchers = [
            SM(r"\s*ensemble\s*\:\s*(?P<ensemble_type>[A-Za-z]+)"),
            SM(r"\s*temperature\s*\:\s*(?P<castep_thermostat_target_temperature>[0-9.]+\.)"),
            SM(r"\s*using\s*\:\s*(?P<castep_barostat_type>[A-Za-z]+\-[A-Za-z]+)\s+"),
            SM(r"\s*with characteristic cell time\s*\:\s*(?P<castep_barostat_tau>[-+0-9.eEd]+)"),
            SM(r"\s*using\s*\:\s*(?P<castep_thermostat_type>[A-Za-z]+\-[A-Za-z]+)\s+"),
            SM(r"\s*with characteristic ionic time\s*\:\s*(?P<castep_thermostat_tau>[-+0-9.eEd]+)"),
            SM(r"\s*time\sstep\s*\:\s*(?P<castep_integrator_dt>[-+0-9.eEdD]+)"),
            SM(r"\s*number of MD steps\s*\:\s*(?P<castep_number_of_steps_requested>[0-9]+)"),
            ])
    ########################################
    # submatcher for section system description
    systemDescriptionSubMatcher = SM(name = "systemDescription",
        startReStr = r"\s*Unit Cell\s*",
        forwardMatch = True,
        sections = ["section_system"],
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
              sections = ["castep_section_atom_positions"],
              subMatchers = [

                 SM(r"\s*a \=\s*(?P<castep_cell_length_a>[\d\.]+)\s*alpha \=\s*(?P<castep_cell_angle_alpha>[\d\.]+)"),
                 SM(r"\s*b \=\s*(?P<castep_cell_length_b>[\d\.]+)\s*beta  \=\s*(?P<castep_cell_angle_beta>[\d\.]+)"),
                 SM(r"\s*c \=\s*(?P<castep_cell_length_c>[\d\.]+)\s*gamma \=\s*(?P<castep_cell_angle_gamma>[\d\.]+)"),
                 SM(r"\s*x\s*(?P<castep_store_atom_labels>[A-Za-z0-9]+\s+[\d\.]+)\s*[0-9]\s*(?P<castep_store_atom_positions>[\d\.]+\s+[\d\.]+\s+[\d\.]+)",
                    endReStr = "\n",
                    repeats = True)

                             ]), # CLOSING castep_section_atom_position

            # atomic positions and cell dimesions
           SM(startReStr = r"\s*Units of ionic velocities is ANG\/PS\s*",
              forwardMatch = True,
              #sections = ["castep_section_atom_position"],
              subMatchers = [
                 SM(r"\s*x\s*[A-Za-z0-9]+\s+[\d\.]+\s*[0-9]\s*(?P<castep_store_atom_ionic_velocities>[-+0-9.eEdD]+\s*[-+0-9.eEdD]+\s*[-+0-9.eEdD]+)",
                    endReStr = "\n",
                    repeats = True)

                             ]), # CLOSING castep_section_atom_positions

                      ]) # CLOSING SM systemDescription

   

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


    singlepointSubMatcher = SM(name = 'single_point',
                startReStr = r"\-*\s*\<\-\-\sSCF\s*",             
                #endReStr = r"BFGS\: finished iteration     0 with enthalpy\= \-2\.14201700E\+002 eV",
                sections = ["section_single_configuration_calculation"],
                forwardMatch = True,
                # endReStr = "\n",
                # endReStr =r"\sStarting\s[A-Z]+\siteration\s*[0-9.]+\s\.\.\.\s*",
                subMatchers = [                     
                    #geomOptimSubMatcher_init, 
                    
                    SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                             endReStr = "\n",
                             repeats = True),
                    SM(sections = ['section_scf_iteration'],
                        startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                        endReStr = "\n",
                        repeats = True),
                    
                    SM(r"Final energy = *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                    SM(r"Final energy\,\s*E\s*= *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                    SM(r"\sTotal energy corrected for finite basis set\s\=\s*(?P<CASTEP_total_energy_corrected_for_finite_basis>[-+0-9.eEdD]+)\s+"),
                    SM(r"Dispersion corrected final energy\*\s=\s*(?P<castep_total_dispersion_corrected_energy__eV>[-+0-9.eEdD]*)"),#total energy including dispersion correction
                    SM(r"Final free energy\s*\(E\-TS\)\s*= *(?P<energy_free__eV>[-+0-9.eEdD]*)"), # matching final converged total free energy
                    SM(r"NB est\. 0K energy\s*\(E\-0\.5TS\)\s*= *(?P<energy_total_T0__eV>[-+0-9.eEdD]*)"), # 0K corrected final SCF energy
                 
                    SM(startReStr = r"\s\*\*\*\*\** Forces \*\*\*\*\**\s*",
                         subMatchers = [
                                    SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                        repeats = True)
                         ]), 

                    SM(#startReStr = r"\sBFGS\:\sFinal\sbulk\smodulus\sunchanged\sfrom\sinitial\svalue\s*",
                        startReStr = r"\s\*\*\*\*\**\sSymmetrised Forces\s\*\*\*\*\**\s*",
                        subMatchers = [
                           SM(r"\s\*\s[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True)
                                      ]),
                     
                    SM(name = 'stresstensor',
                        #startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\** Stress Tensor \*\*\*\*\*\*\*\*\*\*\**\s*",
                        startReStr = r"\s\*\*\*\*\** Stress Tensor \*\*\*\*\**\s*",
                        #GGGrepeats = True,
                        sections = ['castep_section_stress_tensor'],
                        subMatchers = [
                           SM(r"\s*\*\s*[a-z]\s*(?P<castep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True),  
                                      ]), # CLOSING section_stress_tensor
                   
                     
                    SM(name = 'stresstensor',
                        #startReStr = r"\s\*\*\*\*\*\*\*\*\*\ Symmetrised Stress Tensor \*\*\*\*\*\*\*\*\*\*\*\s*",
                        startReStr = r"\s\*\*\*\*\** Symmetrised Stress Tensor \*\*\*\*\**\s*",
                        sections = ['castep_section_stress_tensor'],
                        subMatchers = [
                           SM(r"\s*\*\s*[a-z]\s*(?P<castep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True),  
                                      ]), # CLOSING section_stress_tensor
                    SM(startReStr = r"\s*x\s*MD\sData\:\s*x",
                         subMatchers = [
                            SM(r"\s*x\s*time\s*\:\s*(?P<castep_frame_time_0>[+0-9.eEdD]+)\s*ps\s*x\s*"),
                    ]),

                 ])     
    
    MDSubMatcher = SM(name = 'MD',
                startReStr = r"\sStarting MD iteration\s*[0-9.]+\s\.\.\.\s*",             
                #endReStr = r"BFGS\: finished iteration     0 with enthalpy\= \-2\.14201700E\+002 eV",
                sections = ["castep_section_SCF_iteration_frame"],
                # sections = ["section_single_configuration_calculation"],
                endReStr = r"\s\.\.\.\sfinished MD iteration\s*[0-9.]+\s*",
                repeats =True,
                subMatchers = [                     
                    # SM(r"\s*[0-9]+\s*(?P<castep_SCF_frame_energy>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                    #     endReStr = "\n",
                    #     repeats = True),
                    SM(r"\s*[0-9]+\s*(?P<castep_SCF_frame_energy>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                        endReStr = "\n",
                        repeats = True),                
                    
                    SM(startReStr = r"\s*x\s*MD\sData\:\s*x",
                         subMatchers = [
                            SM(r"\s*x\s*time\s*\:\s*(?P<castep_frame_time>[+0-9.eEdD]+)\s*ps\s*x\s*"),
                    ]),
                    ])
    ########################################
    # Sub matcher for geometry optimisation 
    ########################################

    geomOptimSubMatcher_improving_cycle =  SM (name = 'geometry_optimisation_improving',
            startReStr = r"\s[A-Za-z]+\:\simproving\siteration\s*1\s+",#with\sline\sminimization\s\(lambda\=\s*1\.557176\)\s*",
       
            endReStr = r"\s[A-Za-z]+\:\sfinished iteration\s*1",
       
            subMatchers = [               
              
                        SM(name = 'cellInformation',
                            startReStr = r"\s*Unit Cell\s*",
                            forwardMatch = True,
                            sections = ["castep_section_cell_optim"],
                            subMatchers = [
                                SM(r"\s*(?P<castep_cell_vector_optim>[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s+[-+0-9.eEd]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                                    endReStr = "\n",
                                    repeats = True),

                             ]), # CLOSING castep_section_cell


           # atomic positions and cell dimesions
                        SM(startReStr = r"\s*Lattice parameters",
                            forwardMatch = True,
                            sections = ["castep_section_atom_positions_optim"],
                            subMatchers = [

                                SM(r"\s*a \=\s*(?P<castep_cell_length_a_optim>[\d\.]+)\s*alpha \=\s*(?P<castep_cell_angle_alpha_optim>[\d\.]+)"),
                                SM(r"\s*b \=\s*(?P<castep_cell_length_b_optim>[\d\.]+)\s*beta  \=\s*(?P<castep_cell_angle_beta_optim>[\d\.]+)"),
                                SM(r"\s*c \=\s*(?P<castep_cell_length_c_optim>[\d\.]+)\s*gamma \=\s*(?P<castep_cell_angle_gamma_optim>[\d\.]+)"),

                            ]), # CLOSING castep_section_atom_positions

                        SM(r"\s*x\s*(?P<castep_store_optimised_atom_labels>[A-Za-z]+\s*[0-9]+)\s*(?P<castep_store_optimised_atom_positions>[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)",
                            endReStr = "\n",
                            repeats = True),

                        
                        SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                             endReStr = "\n",
                             repeats = True),

                                   # ]), # CLOSING section_scf_iteration                                      
                        SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                            endReStr = "\n",
                            repeats = True),           
                            
                        SM(r"Final energy = *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                        SM(r"Final energy\,\s*E\s*= *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                 #SM(r"Final free energy\s*\(E\-TS\)\s*= *(?P<castep_energy_free>[-+0-9.eEdD]*)"),
                

                        SM(name = 'stresstensor',
                         #startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\** Stress Tensor \*\*\*\*\*\*\*\*\*\*\**\s*",
                        startReStr = r"\s\*\*\*\*\** Stress Tensor \*\*\*\*\**\s*",
                        forwardMatch = True,
                        #repeats = True,
                        sections = ['castep_section_stress_tensor'],
                        subMatchers = [
                            SM(r"\s*\*\s*[a-z]\s*(?P<castep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                               repeats = True),  
                                       ]), # CLOSING section_stress_tensor

                        
                

                 ])
    geomOptimSubMatcher =  SM (name = 'geometry_optimisation',
            startReStr = r"\sStarting BFGS iteration\s*(?P<CASTEP_geom_iteration_index>[0-9]+)\s\.\.\.\s*",
            sections = ['section_single_configuration_calculation','section_system'],
            endReStr = r"\s*Atomic\sPopulations\s\(Mulliken\)\s*",
            #endReStr = r"\s\[A-Za-z]+\:\sGeometry\soptimization\scompleted\ssuccessfully.\s*",
            repeats = True,
            subMatchers = [               
                        #geomOptimSubMatcher_improving_cycle,
                        SM(name = 'cellInformation',
                            startReStr = r"\s*Unit Cell\s*",
                            forwardMatch = True,
                            sections = ["castep_section_cell_optim"],
                            subMatchers = [
                                SM(r"\s*(?P<castep_cell_vector_optim>[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s+[-+0-9.eEd]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                                    endReStr = "\n",
                                    repeats = True),

                             ]), # CLOSING castep_section_cell


           # atomic positions and cell dimesions
                        SM(startReStr = r"\s*Lattice parameters",
                            forwardMatch = True,
                            sections = ["castep_section_atom_positions_optim"],
                            subMatchers = [

                                SM(r"\s*a \=\s*(?P<castep_cell_length_a_optim>[\d\.]+)\s*alpha \=\s*(?P<castep_cell_angle_alpha_optim>[\d\.]+)"),
                                SM(r"\s*b \=\s*(?P<castep_cell_length_b_optim>[\d\.]+)\s*beta  \=\s*(?P<castep_cell_angle_beta_optim>[\d\.]+)"),
                                SM(r"\s*c \=\s*(?P<castep_cell_length_c_optim>[\d\.]+)\s*gamma \=\s*(?P<castep_cell_angle_gamma_optim>[\d\.]+)"),

                            ]), # CLOSING castep_section_atom_positions

                        SM(r"\s*x\s*(?P<castep_store_optimised_atom_labels>[A-Za-z]+\s*[0-9]+)\s*(?P<castep_store_optimised_atom_positions>[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)",
                            endReStr = "\n",
                            repeats = True),

                        
                        SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                             endReStr = "\n",
                             repeats = True),

                                   # ]), # CLOSING section_scf_iteration                                      
                        SM(sections = ['section_scf_iteration'],
                            startReStr = r"\s*[0-9]+\s*(?P<energy_total_scf_iteration__eV>[-+0-9.eEdD]*)\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[0-9.]*\s*\<\-\-\sSCF\s*",
                            endReStr = "\n",
                            repeats = True),           
                            
                        SM(r"Final energy = *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                        SM(r"Final energy\,\s*E\s*= *(?P<energy_total__eV>[-+0-9.eEdD]*)"), # matching final converged total energy
                        #SM(r"Final free energy\s*\(E\-TS\)\s*= *(?P<energy_free>[-+0-9.eEdD]*)"),
                
                
                        SM(name = 'Forces',
                            startReStr = r"\s\*\*\*\*\** Forces \*\*\*\*\**\s*",
                            subMatchers = [
                                SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                    repeats = True)
                                      ]),

                        
                        SM(name = 'stresstensor',   
                            #startReStr = r"(energy not corrected for finite basis set)\s*",
                            startReStr = r"\s\*\*\*\*\** Stress Tensor \*\*\*\*\**\s*",
                            
                            #forwardMatch = True,
                            sections = ['castep_section_stress_tensor'],
                            subMatchers = [
                                SM(r"\s*\*\s*[a-z]\s*(?P<castep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                    repeats = True),  
                                       ]), # CLOSING section_stress_tensor
                        
                        
                        geomOptimSubMatcher_improving_cycle,
                        SM(r"\s[A-Za-z]+\:\sGeometry\soptimization\scompleted\s(?P<CASTEP_geom_converged>[a-z]+)\.\s*"),
                

                 ])

    geomOptim_finalSubMatcher = SM (name = 'geometry_optimisation_final_configuration',
            startReStr = r"\s[A-Za-z]+\:\sFinal Configuration\:",
            sections = ['section_single_configuration_calculation','section_system'],

            repeats = True,
            subMatchers = [                  
                                
                        SM(name = 'cellInformation',
                                    startReStr = r"\s*Unit Cell\s*",
                                    forwardMatch = True,
                                    sections = ["castep_section_cell_optim"],
                                    subMatchers = [
                                SM(r"\s*(?P<castep_cell_vector_optim>[-+0-9.eEdD]+\s+[-+0-9.eEdD]+\s+[-+0-9.eEd]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                         #SM(r"\s*(?P<castep_cell_vector>[\d\.]+\s+[\d\.]+\s+[\d\.]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                                        endReStr = "\n",
                                    repeats = True),

                                     ]), # CLOSING castep_section_cell


                   # atomic positions and cell dimesions
                        SM(startReStr = r"\s*Lattice parameters",
                            forwardMatch = True,
                            sections = ["castep_section_atom_positions_optim"],
                            subMatchers = [

                                SM(r"\s*a \=\s*(?P<castep_cell_length_a_optim>[\d\.]+)\s*alpha \=\s*(?P<castep_cell_angle_alpha_optim>[\d\.]+)"),
                                SM(r"\s*b \=\s*(?P<castep_cell_length_b_optim>[\d\.]+)\s*beta  \=\s*(?P<castep_cell_angle_beta_optim>[\d\.]+)"),
                                SM(r"\s*c \=\s*(?P<castep_cell_length_c_optim>[\d\.]+)\s*gamma \=\s*(?P<castep_cell_angle_gamma_optim>[\d\.]+)"),

                                     ]), # CLOSING castep_section_atom_positions

                                SM(r"\s*x\s*(?P<castep_store_optimised_atom_labels>[A-Za-z]+\s*[0-9]+)\s*(?P<castep_store_optimised_atom_positions>[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)",
                                    endReStr = "\n",
                                    repeats = True),
                            

                        
                        SM(r"\s[A-Za-z]+\:\sFinal\sEnthalpy\s*\=\s(?P<CASTEP_enthalpy>[-+0-9.eEdD]+)"),        
                        
                        SM(#startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\* Forces \*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\s*",
                                startReStr = r"\s\*\*\*\*\** Forces \*\*\*\*\**\s*",
                                
                                subMatchers = [
                                   SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                      repeats = True)
                                              ]),
                        SM(#startReStr = r"\sBFGS\:\sFinal\sbulk\smodulus\sunchanged\sfrom\sinitial\svalue\s*",
                                startReStr = r"\s\*\*\*\*\**\sSymmetrised\sForces\s\*\*\*\*\**\s*",
                                repeats = True,
                                subMatchers = [
                                   SM(r"\s\*\s[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                      repeats = True)
                                              ]), 

                        
                        SM(name = 'stresstensor',
                                #startReStr = r"\s\*\*\*\*\*\*\*\*\*\ Symmetrised Stress Tensor \*\*\*\*\*\*\*\*\*\*\*\s*",
                                startReStr = r"\s\*\*\*\*\** Symmetrised Stress Tensor \*\*\*\*\**\s*",
                                sections = ['castep_section_stress_tensor'],
                                subMatchers = [
                                   SM(r"\s*\*\s*[a-z]\s*(?P<castep_store_stress_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                      repeats = True),  
                                              ]), # CLOSING section_stress_tensor

               

                          ])
                
    

    Mulliken_SubMatcher = SM (name = 'Mulliken population analysis',
            startReStr = r"\s*Atomic Populations\s\(Mulliken\)\s*",
            sections = ['castep_section_population_analysis'],
            #endReStr = r"\s*Bond\s*Population\s*Length\s\(A\)\s*",
            repeats = True,
            subMatchers = [ 
                SM(r"\s*[a-zA-Z]+\s*[0-9.]+\s*(?P<castep_orbital_contributions>[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+\s*[-\d\.]+)\s*[-+0-9.eEdD]+\s*(?P<castep_mulliken_charge_store>[-+0-9.eEdD]+)\s*",       
                    endReStr = r"\s*Bond\s*Population\s*Length\s\(A\)\s*",
                    repeats = True),
                 ]) 
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

               
                basisSetCellAssociatedSubMatcher, # section_basis_set_cell_dependent
                
                ElectronicParameterSubMatcher,
                
                MDParameterSubMatcher,
                
                GeomOptimParameterSubMatcher,
                
                phononCalculationSubMatcher,
                
                systemDescriptionSubMatcher, # section_system subMatcher
        
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
               

                SM(r"\s*Point group of crystal\s\=\s*[0-9.]+\:\s(?P<castep_crystal_point_group>[a-zA-Z0-9.]+)"),
                SM(r"\s*Space group of crystal\s\=\s*[0-9.]+\:\s(?P<castep_space_group>[a-zA-Z0-9.]+)"),

          ############ CASTEP-specific van der Waals method parameters #############################     
                
                SM(name = "van der Waals castep TS",
                   startReStr = r"\s*Dispersion\-correction scheme\s\:\s+",
                   forwardMatch = True,
                   sections = ["castep_section_van_der_Waals_parameters"],
                   subMatchers = [
                        ######## Method TS #######
                        SM(r"\s*Parameter sR\s*\: *(?P<Parameter_sR> [0-9.]+)"),
                        SM(r"\s*Parameter d\s*\: *(?P<Parameter_d> [0-9.]+)"),
                        ######## Method OBS #######
                        SM(r"\s*Parameter lambda\s*\: *(?P<Parameter_LAMBDA> [0-9.]+)"),
                        SM(r"\s*Parameter n\s*\: *(?P<Parameter_n> [0-9.]+)"),
               
                        ######## Method G06 #######
                        SM(r"\s*Parameter s6\s*\: *(?P<Parameter_s6> [0-9.]+)"),
                        SM(r"\s*Parameter d\s*\: *(?P<Parameter_d> [0-9.]+)")
                              ]), # CLOSING van der waals castep parameters
              
                #geomOptimSubMatcher_init,
             
                
                                     
                SM(r"Calculating total energy with cut\-off of  (?P<castep_basis_set_planewave_cutoff_iteration_0>[0-9.]+)",
                               repeats = True,),
                
                singlepointSubMatcher,
                
                MDSubMatcher,

                SM(name = "Vibrational_frequencies",
                    sections = ["castep_section_vibrational_frequencies"],
                    startReStr = r"\s\+\s*q\-pt\=\s*[0-9.]+",
                    #startReStr = r"\s\+\s*Vibrational Frequencies\s*\+\s*",
                    repeats = True,
                    subMatchers = [
                        SM(r"\s\+\s*(?P<castep_n_iterations_phonons>[0-9.]+)\s*(?P<castep_vibrationl_frequencies_store>[-\d\.]+)\s*(?P<castep_ir_store>[a-z])\s*(?P<castep_ir_intensity_store>[-\d\.]+)\s*[A-Z]\s*(?P<castep_raman_activity_store>[-\d\.]+)\s*[A-Z]\s*\+\s*",       
                            repeats = True,
                            endReStr = r"\s\+\s\.*\s\+\s*",
                             ),
                        SM(r"\s\+\s*(?P<castep_n_iterations_phonons>[0-9.]+)\s*(?P<castep_vibrationl_frequencies_store>[-\d\.]+)\s*(?P<castep_ir_store>[a-z])\s*(?P<castep_ir_intensity_store>[-\d\.]+)\s*[A-Z]\s*[A-Z]\s*\+\s*",       
                            repeats = True,
                            endReStr = r"\s\+\s\.*\s\+\s*",
                             ),
                        SM(r"\s\+\s*(?P<castep_n_iterations_phonons>[0-9.]+)\s*(?P<castep_vibrationl_frequencies_store>[-\d\.]+)\s*(?P<castep_ir_store>[a-z])\s*\+\s*",          
                            repeats = True,
                            endReStr = r"\s\+\s\.*\s\+\s*",
                            ),
                          ]),
                
                bandStructureSubMatcher,  # band structure subMatcher
                                       
                SM(name = 'ramantensor',
                        #startReStr = r"\s\*\*\*\*\*\*\*\*\*\*\** Stress Tensor \*\*\*\*\*\*\*\*\*\*\**\s*",
                        startReStr = r"\s\+\sMode number\:\s*[0-9.]+\sRaman\stensor\s*Depolarisation Ratio\s*\+\s*",
                        endReStr = r"\s\+\s*\+\s*",
                        repeats = True,
                        sections = ['section_castep_raman_tensor'],
                        subMatchers = [
                           SM(r"\s\+\s*(?P<castep_store_raman_tensor>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True),  
                                      ]), # CLOSIN

                # SM(startReStr = r"\s\*\*\*\*\** Forces \*\*\*\*\**\s*",
                #     subMatchers = [
                #            SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                #               repeats = True)
                #                       ]),

                     
                SM(startReStr = r"\s\*\*\*\*\**\sSymmetrised Forces\s\*\*\*\*\**\s*",
                        subMatchers = [
                           SM(r"\s\*\s[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces_band>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                              repeats = True)
                                      ]),
                    
                geomOptimSubMatcher,
                
                geomOptim_finalSubMatcher,            
                
                Mulliken_SubMatcher,
                
                SM(name = 'calc_time',
                    startReStr = r" A BibTeX formatted list of references used in this run has been written to",
                    sections = ['castep_section_time'],
                    subMatchers = [
                        SM(r"Initialisation time =\s*(?P<castep_initialisation_time>[0-9.]*)"), # matching final converged total energy
                        SM(r"Calculation time    =\s*(?P<castep_calculation_time>[0-9.]*)"),
                        SM(r"Finalisation time   =\s*(?P<castep_finalisation_time>[0-9.]*)"),
                                      ]), # CLOSING section_castep_time
                                        
                  
                
                         ])       
                
 
        ]) # CLOSING SM NewRun 



def get_cachingLevelForMetaName(metaInfoEnv):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = { 'castep_ir_intensity_store' : CachingLevel.Cache,
                                'castep_ir_store' : CachingLevel.Cache,
                                'castep_vibrationl_frequencies_store' : CachingLevel.Cache,
                                'castep_n_iterations_phonons' : CachingLevel.Cache,
                                'castep_mulliken_charge_store' : CachingLevel.Cache,
                                'castep_orbital_contributions' : CachingLevel.Cache,
                                'castep_section_cell_optim': CachingLevel.Cache,
                                'castep_section_atom_positions_optim' : CachingLevel.Cache,
                                'castep_section_eigenvalues':CachingLevel.Cache,
                                'castep_section_k_points':CachingLevel.Cache,
                                'castep_section_k_band':CachingLevel.Cache,
                                'band_energies' : CachingLevel.Cache,
                                #'band_k_points' : CachingLevel.Cache,
                                'castep_basis_set_planewave_cutoff' : CachingLevel.Cache,
                                'eigenvalues_values': CachingLevel.Cache,
                                'eigenvalues_kpoints':CachingLevel.Cache,
                                'castep_number_of_electrons_store':CachingLevel.Cache,
                                 'castep_frame_time':CachingLevel.Cache,
                                'castep_section_SCF_iteration_frame':CachingLevel.Cache,
                                'castep_SCF_frame_energy':CachingLevel.Cache}

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








