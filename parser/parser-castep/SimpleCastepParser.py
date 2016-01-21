import setup_paths
import numpy as np
import math
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.local_meta_info import loadJsonFile, InfoKindEl
from nomadcore.caching_backend import CachingLevel
import re, os, sys, json, logging



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
        self.at_nr                             = 0
        self.atom_label                        = []
        self.atom_forces                       = []
        self.castep_atom_position              = []
        self.atom_position                     = []
        self.a                                 = []
        self.b                                 = []
        self.c                                 = []
        self.alpha                             = []
        self.beta                              = []
        self.gamma                             = []
        self.volume                            = 0

        self.energy_total_scf_iteration_list   = []
        self.scfIterNr                         = []
        self.ecut                              = []
        self.k_count                           = 0
        self.k_nr                              = 0
        self.e_nr                              = 0
        self.eigenvalues_kpoints               = []
        self.eigenvalues_eigenvalues           = []



# Translating the XC functional name to the NOMAD standard (as seen in the CP2K parser, thanks!)
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
            " Perdew Burke Ernzerhof": "GGA_C_PBE",
            " Local Density Approximation": "LDA_C_PZ",
            " Perdew Wang (1991)": "GGA_X_PW91",
            " revised Perdew Burke Ernzerhof": "GGA_X_RPBE",
            " PBE with Wu-Cohen exchange": "GGA_X_WC",
            " PBE for solids (2008)": "GGA_X_PBE_SOL",
        }

        # Define a mapping for the relativistic treatments
        relativistic_map = {
            " Koelling-Harmon": "scalar_relativistic"
        }

        # Match each castep functional name and sort the matches into a list
        functionals = []

        for name in functional_names:
            match = functional_map.get(name)
            if match:
                functionals.append(match)
        functionals = "_".join(sorted(functionals))

        # Push the functional string into the backend
        backend.addValue('XC_functional', functionals)

        # Match each castep relativity treatment name and sort the matches into a list
        relativistic = []

        for name in relativistic_names:
            match = relativistic_map.get(name)
            if match:
                relativistic.append(match)
        relativistic = "_".join(sorted(relativistic))

        # Push the relativistic treatment string into the backend
        backend.addValue('relativity_method', relativistic)



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
        backend.addArrayValues('simulation_cell', np.asarray(self.cell), unit='angstrom')



# Storing the total energy of each SCF iteration in an array
    def onClose_section_scf_iteration(self, backend, gIndex, section):
        """trigger called when _section_scf_iteration is closed"""
        # get cached values for energy_total_scf_iteration
        ev = section['energy_total_scf_iteration']
        self.scfIterNr = len(ev)
        self.energy_total_scf_iteration_list.append(ev)

        backend.addArrayValues('energy_total_scf_iteration_list', np.asarray(self.energy_total_scf_iteration_list))
        backend.addValue('scf_dft_number_of_iterations', self.scfIterNr)



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



# Here we add basis set name and kind for the plane wave code
    def onClose_section_basis_set_cell_associated(self, backend, gIndex, section):
        ecut_str = section['castep_basis_set_plan_wave_cutoff']
        self.ecut = float(ecut_str[0])
        eVtoRy = 0.073498618
        ecut_str_name = str(int(round(eVtoRy*self.ecut)))

        basis_set_kind = 'plane_waves'
        basis_set_name = 'PW_'+ecut_str_name
        backend.addValue('basis_set_plan_wave_cutoff', self.ecut)
        backend.addValue('basis_set_cell_associated_kind', basis_set_kind)
        backend.addValue('basis_set_cell_associated_name', basis_set_name)



######################################################################################
################ Triggers on closure section_system_description ######################
######################################################################################

    def onClose_section_system_description(self, backend, gIndex, section):
        """trigger called when _section_system_description is closed"""


# Processing forces acting on atoms (final converged forces)
        #get cached values of castep_store_atom_forces
        f_st = section['castep_store_atom_forces']
        self.at_nr = len(f_st)
        for i in range(0, self.at_nr):
            f_st[i] = f_st[i].split()
            f_st[i] = [float(j) for j in f_st[i]]
            f_st_int = f_st[i]
            self.atom_forces.append(f_st_int)
        backend.addArrayValues('atom_forces', np.asarray(self.atom_forces))



# Processing the atom positions in fractionary coordinates (as given in the CASTEP output)
        #get cached values of castep_store_atom_position
        pos = section['castep_store_atom_position']
        for i in range(0, self.at_nr):
            pos[i] = pos[i].split()
            pos[i] = [float(j) for j in pos[i]]
            self.castep_atom_position.append(pos[i])
        backend.addArrayValues('castep_atom_position', np.asarray(self.castep_atom_position))



# Processing the atom labels
        #get cached values of castep_store_atom_label
        lab = section['castep_store_atom_label']
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
        backend.addArrayValues('atom_position', np.asarray(self.atom_position))



# Backend add the total number of atoms in the simulation cell
        backend.addValue('number_of_atoms', self.at_nr)



######################################################################################
###################### Storing k points and band energies ############################
######################################################################################

# Storing the k point coordinates
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
            self.eigenvalues_kpoints.append(k_st_int)



# Storing the eigenvalues
    def onClose_castep_section_eigenvalues(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""

        #get cached values of castep_store_k_points
        e_st = section['castep_store_eigenvalues']
        self.e_nr = len(e_st)
        self.eigenvalues_eigenvalues.append(e_st)



# Keeping only the last arrays update
    def onClose_section_eigenvalues(self, backend, gIndex, section):
        """trigger called when _section_eigenvalues"""
        backend.addArrayValues('eigenvalues_kpoints', np.asarray(self.eigenvalues_kpoints))
        backend.addArrayValues('eigenvalues_eigenvalues', np.asarray(self.eigenvalues_eigenvalues))
# Backend add the number of k points and eigenvalues
        backend.addValue('number_of_eigenvalues_kpoints', self.k_nr)
        backend.addValue('number_of_eigenvalues', self.e_nr)



################################################################################################################################################################
################################################################################################################################################################
################################################################################################################################################################
######################  MAIN PARSER STARTS HERE  ###############################################################################################################
################################################################################################################################################################
############################################################################ CASTEP.Parser Version 1.0 #########################################################
################################################################################################################################################################

# main Parser
mainFileDescription = SM(name = 'root',
                         weak = True,
                         startReStr = "",
                         subMatchers = [

                            SM(name = 'newRun',
                            startReStr = r"\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*",
                            #repeats = True,
                            required = True,
                            forwardMatch = True,
                            sections   = ['section_run'],
                            subMatchers = [

                                SM(name = 'ProgramHeader',
                                   startReStr = r"\s\|\s*CCC\s*AA\s*SSS\s*TTTTT\s*EEEEE\s*PPPP\s*\|\s*",
                                   subMatchers = [

                                    SM(r"\s\|\sWelcome to Academic Release\s(?P<program_name>[a-zA-Z]+)* version *(?P<program_version>[0-9a-zA-Z_.]*)"),
                                    SM(r"\sCompiled for *(?P<program_compilation_host>[-a-zA-Z0-9._]*)\son\s(?P<castep_program_compilation_date>[a-zA-Z,\s0-9]*)\s *(?P<castep_program_compilation_time>[0-9:]*)"),
                                    #SM(r"\sRun started\:\s[A-Za-z,]\s*(?P<castep_program_execution_date> [\d\]+\s[A-Za-z]+\s[\d\]+)\s *(?P<castep_program_execution_time>[0-9:]*)"),
                                    SM(r"\sCompiler\: *(?P<castep_compiler>[a-zA-Z\s0-9.]*)"),
                                    SM(r"\sMATHLIBS\: *(?P<castep_maths_library>[a-zA-Z0-9.() ]*)\s*"),
                                    SM(r"\sFFT Lib \: *(?P<castep_fft_library>[a-zA-Z0-9.() ]*)\s*"),
                                    SM(r"\sFundamental constants values\: *(?P<castep_constants_reference>[a-zA-Z0-9.() ]*)\s*"),

                                  ]), # CLOSING SM ProgramHeader


                                SM(startReStr = r"\sCalculation not parallelised\.",
                                   sections = ["section_system_description"],
                                   subMatchers = [

                                       # XC functional
                                       SM(name = 'XCMethods',
                                          startReStr = r"\susing functional\s*\:",
                                          forwardMatch = True,
                                          sections = ["section_method"],
                                          subMatchers = [

                                              SM(name = "castepXC",
                                                 startReStr = r"\susing functional\s*\:",
                                                 forwardMatch = True,
                                                 sections = ["castep_section_functionals"],
                                                 otherMetaInfo = ["XC_functional"],
                                                 dependencies = {"XC_functional": ["castep_functional_name"]},
                                                 subMatchers = [

                                                    SM(r"\susing functional\s*\: *(?P<castep_functional_name> [A-Za-z0-9() ]*)"),
                                                    SM(r"\srelativistic treatment\s*\: *(?P<castep_relativity_treatment_scf> [A-Za-z0-9() -]*)")

                                                 ]), # CLOSING castep_section_functionals


                                          ]), # CLOSING section_method


                                       # cell information
                                       SM(name = 'planeWave basis set',
                                          startReStr = r"\sbasis set accuracy\s*",
                                          forwardMatch = True,
                                          sections = ["section_basis_set_cell_associated"],
                                          subMatchers = [

                                              SM(r"\splane wave basis set cut\-off\s*\:\s*(?P<castep_basis_set_plan_wave_cutoff>[0-9.]+)")

                                          ]), # CLOSING section_basis_set_cell_associated


                                       # cell information
                                       SM(name = 'cellInformation',
                                          startReStr = r"\s*Unit Cell\s*",
                                          forwardMatch = True,
                                          sections = ["castep_section_cell"],
                                          subMatchers = [

                                            SM(r"\s*(?P<castep_cell_vector>[\d\.]+\s+[\d\.]+\s+[\d\.]+) \s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*\s*[-+0-9.eEdD]*",
                                               endReStr = "\n",
                                               repeats = True),

                                          ]), # CLOSING castep_section_cell


                                       # atomic positions and cell dimesions
                                       SM(startReStr = r"\s*Lattice parameters",
                                          sections = ["castep_section_atom_position"],
                                          subMatchers = [

                                            SM(r"\s*a \=\s*(?P<castep_cell_length_a>[\d\.]+)\s*alpha \=\s*(?P<castep_cell_angle_alpha>[\d\.]+)"),
                                            SM(r"\s*b \=\s*(?P<castep_cell_length_b>[\d\.]+)\s*beta  \=\s*(?P<castep_cell_angle_beta>[\d\.]+)"),
                                            SM(r"\s*c \=\s*(?P<castep_cell_length_c>[\d\.]+)\s*gamma \=\s*(?P<castep_cell_angle_gamma>[\d\.]+)"),
                                            SM(r"\s*x\s*(?P<castep_store_atom_label>[A-Za-z0-9]+\s+[\d\.]+)\s*[0-9]\s*(?P<castep_store_atom_position>[\d\.]+\s+[\d\.]+\s+[\d\.]+)",
                                               endReStr = "\n",
                                               repeats = True)

                                          ]), # CLOSING castep_section_atom_position


                                       # SCF Single Point Evaluation (energies and forces)
                                       SM(name = 'SinglePointEvaluation',
                                          startReStr = r"SCF\sloop\s*Energy\s*Energy\sgain\s*Timer\s*<\-\-\sSCF\s*",
                                          #repeats = True,
                                          forwardMatch = True,
                                          sections = ['section_single_configuration_calculation'],
                                          subMatchers = [

                                              SM(name = 'ScfIterations',
                                                 startReStr = r"SCF\sloop\s*Energy\s*Energy\sgain\s*Timer\s*<\-\-\sSCF\s*",
                                                 sections = ['section_scf_iteration'],
                                                 subMatchers = [

                                                    SM(r"\s*[1-9]\s*(?P<energy_total_scf_iteration>[-+0-9.eEdD]*)", repeats = True),

                                                 ]), # CLOSING section_scf_iteration


                                          SM(r"Final energy = *(?P<energy_total>[-+0-9.eEdD]*)"), # macthing final coverged total energy

                                          #  Band Structure Calculation: here we match and parse k points and eigenvalues
                                          SM(startReStr = "\s*\+\s*B A N D   S T R U C T U R E   C A L C U L A T I O N\s*",
                                                 forwardMatch = True,
                                                 sections = ["section_eigenvalues_group"],
                                                 subMatchers = [

                                                     SM(startReStr = "\s*\+\s*B A N D   S T R U C T U R E   C A L C U L A T I O N\s*",
                                                        sections = ["section_eigenvalues"],
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

                                                                    SM(name = 'Eigen',
                                                                       startReStr = r"\s*\+\s*\+\s*",
                                                                       sections = ['castep_section_eigenvalues'],
                                                                       repeats = True,
                                                                       subMatchers = [
                                                                          # Matching eigenvalues
                                                                          SM(r"\s*\+\s*[0-9]+\s*(?P<castep_store_eigenvalues>\s+[-\d\.]+)",
                                                                             repeats = True)

                                                                       ]), # CLOSING castep_section_eigenvalues


                                                            ]), # CLOSING castep_section_k_points


                                                        ]), # CLOSING section_eigenvalues


                                                ]), # CLOSING section_eigenvalues_group


                                          SM(r"\s*\*\s*[A-Za-z]+\s*[0-9]\s*(?P<castep_store_atom_forces>[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+)",
                                             endReStr = "\n",
                                             repeats = True)


                                          ]), # CLOSING SM SinglePointEvaluation


                                ]), # CLOSING section_system_description


    ]), # CLOSING SM newRun
])






# loading metadata from nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json
metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/castep.nomadmetainfo.json"))
metaInfoEnv, warnings = loadJsonFile(filePath = metaInfoPath, dependencyLoader = None, extraArgsHandling = InfoKindEl.ADD_EXTRA_ARGS, uri = None)
# dictionary for functions to call on close of a section metaname
onClose = {}
# parser info
parserInfo = {'name':'castep-parser', 'version': '1.0'}
# adjust caching of metadata
cachingLevelForMetaName = {'energy_total': CachingLevel.Cache,
                           'energy_total_scf_iteration_list': CachingLevel.Forward,
                           'atom_forces': CachingLevel.Forward,
                           'atom_label': CachingLevel.Forward,
                           'number_of_atoms': CachingLevel.Forward,
                           'castep_atom_position': CachingLevel.Forward,
                           'atom_position': CachingLevel.Forward,
                           'basis_set_plan_wave_cutoff': CachingLevel.Forward,
                           'castep_basis_set_plan_wave_cutoff': CachingLevel.Cache,
                           'basis_set_cell_associated_kind': CachingLevel.Forward,
                           'basis_set_cell_associated_name': CachingLevel.Forward,

                           'eigenvalues_kpoints': CachingLevel.Forward,
                           'eigenvalues_eigenvalues': CachingLevel.Forward,
                           'number_of_eigenvalues_kpoints': CachingLevel.Forward,
                           'number_of_eigenvalues': CachingLevel.Forward,
                           'castep_store_k_points': CachingLevel.Cache,
                           'castep_store_eigenvalues': CachingLevel.Cache,
                           'castep_store_atom_label': CachingLevel.Cache,
                           }

if __name__ == "__main__":
    mainFunction(mainFileDescription, metaInfoEnv, parserInfo, superContext = CastepParserContext(), cachingLevelForMetaName = cachingLevelForMetaName, onClose = onClose,
                 defaultSectionCachingLevel = False)
