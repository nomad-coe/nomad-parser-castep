#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import numpy as np            # pylint: disable=unused-import
import typing                 # pylint: disable=unused-import
from nomad.metainfo import (  # pylint: disable=unused-import
    MSection, MCategory, Category, Package, Quantity, Section, SubSection, SectionProxy,
    Reference
)
from nomad.metainfo.legacy import LegacyDefinition

from nomad.datamodel.metainfo import public
from nomad.datamodel.metainfo import common

m_package = Package(
    name='castep_nomadmetainfo_json',
    description='None',
    a_legacy=LegacyDefinition(name='castep.nomadmetainfo.json'))


class x_castep_section_vibrational_frequencies(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_vibrational_frequencies'))

    x_castep_number_vibrational_frequencies = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of vibration frequenices
        ''',
        a_legacy=LegacyDefinition(name='x_castep_number_vibrational_frequencies'))

    x_castep_vibrational_frequencies = Quantity(
        type=np.dtype(np.float64),
        shape=['x_castep_number_vibrational_frequencies'],
        description='''
        Vibration Frequenices (cm-1)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_vibrational_frequencies'))

    x_castep_vibrationl_frequencies_store = Quantity(
        type=str,
        shape=[],
        description='''
        Vibration Frequenices (cm-1)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_vibrationl_frequencies_store'))

    x_castep_ir_store = Quantity(
        type=str,
        shape=[],
        description='''
        Irreducible representation in the Point Group
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ir_store'))

    x_castep_ir = Quantity(
        type=str,
        shape=['x_castep_number_vibrational_frequencies'],
        description='''
        Irreducible representation in the Point Group
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ir'))

    x_castep_raman_activity = Quantity(
        type=np.dtype(np.float64),
        shape=['x_castep_number_vibrational_frequencies'],
        description='''
        Raman activity (A**4/amu)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_raman_activity'))

    x_castep_raman_active = Quantity(
        type=str,
        shape=['x_castep_number_vibrational_frequencies'],
        description='''
        Raman activity (A**4/amu)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_raman_active'))

    x_castep_n_iterations_phonons = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of iterations in phonons
        ''',
        a_legacy=LegacyDefinition(name='x_castep_n_iterations_phonons'))

    x_castep_raman_activity_store = Quantity(
        type=str,
        shape=[],
        description='''
        Raman activity (A**4/amu)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_raman_activity_store'))

    x_castep_ir_intensity = Quantity(
        type=np.dtype(np.float64),
        shape=['x_castep_number_vibrational_frequencies'],
        description='''
        IR intensities (D/A)**2/amu
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ir_intensity'))

    x_castep_ir_intensity_store = Quantity(
        type=str,
        shape=[],
        description='''
        IR intensities (D/A)**2/amu
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ir_intensity_store'))


class x_castep_section_band_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_band_parameters'))

    x_castep_band_n_bands = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        x_castep_band_n_bands
        ''',
        a_legacy=LegacyDefinition(name='x_castep_band_n_bands'))

    x_castep_band_conv_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x_castep_band_conv_tolerance
        ''',
        a_legacy=LegacyDefinition(name='x_castep_band_conv_tolerance'))

    x_castep_band_n_iterations = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        x_castep_band_n_iterations
        ''',
        a_legacy=LegacyDefinition(name='x_castep_band_n_iterations'))

    x_castep_band_max_cg = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x_castep_band_max_cg
        ''',
        a_legacy=LegacyDefinition(name='x_castep_band_max_cg'))


class x_castep_section_core_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_core_parameters'))

    x_castep_core_spectra_n_bands = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        x_castep_core_spectra_n_bands
        ''',
        a_legacy=LegacyDefinition(name='x_castep_core_spectra_n_bands'))

    x_castep_core_spectra_conv_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x_castep_core_spectra_conv_tolerance
        ''',
        a_legacy=LegacyDefinition(name='x_castep_core_spectra_conv_tolerance'))


class x_castep_section_ts_scf_iteration(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_scf_iteration'))

    x_castep_scf_ts_iteration_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        SCF_ts_energy
        ''',
        a_legacy=LegacyDefinition(name='x_castep_scf_ts_iteration_energy'))

    x_castep_scf_ts_iteration_energy_change = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        SCF_ts_energy change
        ''',
        a_legacy=LegacyDefinition(name='x_castep_scf_ts_iteration_energy_change'))

    x_castep_scf_ts_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        SCF_ts_energy time
        ''',
        a_legacy=LegacyDefinition(name='x_castep_scf_ts_time'))


class x_castep_section_tddft(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_tddft'))

    x_castep_tddft_iteration = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Iteration number
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_iteration'))

    x_castep_wall_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Iteration wall time (s)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_wall_time'))

    x_castep_state_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        state number
        ''',
        a_legacy=LegacyDefinition(name='x_castep_state_number'))

    x_castep_state_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        state energy
        ''',
        a_legacy=LegacyDefinition(name='x_castep_state_energy'))

    x_castep_state_energy_error = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        state energy error
        ''',
        a_legacy=LegacyDefinition(name='x_castep_state_energy_error'))

    x_castep_tddft_calculation_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        calculation time
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_calculation_time'))


class x_castep_section_atom_ionic_velocities(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_atom_ionic_velocities'))


class x_castep_section_atom_positions_optim(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_atom_positions_optim'))

    x_castep_cell_angle_alpha_optim = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Simulation cell angle alpha
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_angle_alpha_optim'))

    x_castep_cell_angle_beta_optim = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Simulation cell angle beta
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_angle_beta_optim'))

    x_castep_cell_angle_gamma_optim = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Simulation cell angle gamma
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_angle_gamma_optim'))

    x_castep_cell_length_a_optim = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        a unit cell edge length
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_length_a_optim'))

    x_castep_cell_length_b_optim = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        b unit cell edge length
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_length_b_optim'))

    x_castep_cell_length_c_optim = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        c unit cell edge length
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_length_c_optim'))


class x_castep_section_atom_positions(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_atom_positions'))

    x_castep_cell_angle_alpha = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Simulation cell angle alpha
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_angle_alpha'))

    x_castep_cell_angle_beta = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Simulation cell angle beta
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_angle_beta'))

    x_castep_cell_angle_gamma = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Simulation cell angle gamma
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_angle_gamma'))

    x_castep_cell_length_a = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        a unit cell edge length
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_length_a'))

    x_castep_cell_length_b = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        b unit cell edge length
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_length_b'))

    x_castep_cell_length_c = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        c unit cell edge length
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_length_c'))


class x_castep_section_cell_optim(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_cell_optim'))

    x_castep_cell_vector_optim = Quantity(
        type=str,
        shape=[],
        description='''
        Temporay storage for cell vectors
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_vector_optim'))


class x_castep_section_cell(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_cell'))

    x_castep_cell_vector = Quantity(
        type=str,
        shape=[],
        description='''
        Temporay storage for cell vectors
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_vector'))


class x_castep_section_collect_scf_eigenvalues(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_collect_scf_eigenvalues'))


class x_castep_section_eigenvalues_1(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_eigenvalues_1'))

    x_castep_store_eigenvalues_1 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storing eigenvalues
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_eigenvalues_1'))


class x_castep_section_eigenvalues(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_eigenvalues'))

    x_castep_store_eigenvalues = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Temporary storing eigenvalues
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_eigenvalues'))


class x_castep_section_functional_definition(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_functional_definition'))

    x_castep_functional_type = Quantity(
        type=str,
        shape=[],
        description='''
        XC functional definition in CASTEP convention
        ''',
        a_legacy=LegacyDefinition(name='x_castep_functional_type'))

    x_castep_functional_weight = Quantity(
        type=str,
        shape=[],
        description='''
        XC functional weight in CASTEP convention
        ''',
        a_legacy=LegacyDefinition(name='x_castep_functional_weight'))


class x_castep_section_functionals(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_functionals'))

    x_castep_functional_name = Quantity(
        type=str,
        shape=[],
        description='''
        XC functional in CASTEP convention
        ''',
        a_legacy=LegacyDefinition(name='x_castep_functional_name'))


class x_castep_section_population_analysis(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_population_analysis'))

    x_castep_mulliken_atom_index = Quantity(
        type=str,
        shape=[],
        description='''
        Mulliken_atom_index
        ''',
        a_legacy=LegacyDefinition(name='x_castep_mulliken_atom_index'))

    x_castep_mulliken_atom = Quantity(
        type=str,
        shape=[],
        description='''
        Mulliken_atom kind
        ''',
        a_legacy=LegacyDefinition(name='x_castep_mulliken_atom'))

    x_castep_orbital_contributions = Quantity(
        type=str,
        shape=[],
        description='''
        Mulliken_contributions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_orbital_contributions'))

    x_castep_orbital_s = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Mulliken_contribution_orbital s
        ''',
        a_legacy=LegacyDefinition(name='x_castep_orbital_s'))

    x_castep_orbital_p = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Mulliken_contribution_orbital p
        ''',
        a_legacy=LegacyDefinition(name='x_castep_orbital_p'))

    x_castep_orbital_d = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Mulliken_contribution_orbital d
        ''',
        a_legacy=LegacyDefinition(name='x_castep_orbital_d'))

    x_castep_orbital_f = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Mulliken_contribution_orbital f
        ''',
        a_legacy=LegacyDefinition(name='x_castep_orbital_f'))

    x_castep_total_orbital = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Mulliken_total_contribution
        ''',
        a_legacy=LegacyDefinition(name='x_castep_total_orbital'))

    x_castep_mulliken_charge_store = Quantity(
        type=str,
        shape=[],
        description='''
        Mulliken_charges
        ''',
        a_legacy=LegacyDefinition(name='x_castep_mulliken_charge_store'))

    x_castep_mulliken_charge = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms'],
        description='''
        Mulliken_charges
        ''',
        a_legacy=LegacyDefinition(name='x_castep_mulliken_charge'))


class x_castep_section_geom_optimisation_method(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_geom_optimisation_method'))

    x_castep_geometry_optim_method = Quantity(
        type=str,
        shape=[],
        description='''
        Determines optimisation method used
        ''',
        a_legacy=LegacyDefinition(name='x_castep_geometry_optim_method'))


class x_castep_section_optics_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_optics_parameters'))

    x_castep_optics_n_bands = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        optics_number_of_bands
        ''',
        a_legacy=LegacyDefinition(name='x_castep_optics_n_bands'))

    x_castep_optics_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        optics_band_convergence_tolerance
        ''',
        a_legacy=LegacyDefinition(name='x_castep_optics_tolerance'))


class x_castep_section_tddft_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_tddft_parameters'))

    x_castep_tddft_n_excited_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number of excited states
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_n_excited_states'))

    x_castep_tddft_n_states_forces = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number of states for forces
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_n_states_forces'))

    x_castep_tddft_state_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        tolerance (eV)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_state_tolerance'))

    x_castep_tddft_state_tolerance_window = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        tolerance window iterations
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_state_tolerance_window'))

    x_castep_tddft_max_iter = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        max number iterations
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_max_iter'))

    x_castep_tddft_extra_states = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        number of extra states
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_extra_states'))

    x_castep_tddft_functional = Quantity(
        type=str,
        shape=[],
        description='''
        tddft functional
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_functional'))

    x_castep_tddft_method = Quantity(
        type=str,
        shape=[],
        description='''
        tddft method
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_method'))

    x_castep_tddft_eigenmethod = Quantity(
        type=str,
        shape=[],
        description='''
        tddft eigenmethod
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_eigenmethod'))

    x_castep_tddft_approximation = Quantity(
        type=str,
        shape=[],
        description='''
        tddft approximation
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_approximation'))

    x_castep_tddft_position_op = Quantity(
        type=str,
        shape=[],
        description='''
        tddft position operator
        ''',
        a_legacy=LegacyDefinition(name='x_castep_tddft_position_op'))


class x_castep_section_phonons(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_phonons'))

    x_castep_phonon_method = Quantity(
        type=str,
        shape=[],
        description='''
        Phonon calculation method
        ''',
        a_legacy=LegacyDefinition(name='x_castep_phonon_method'))

    x_castep_DFPT_solver_method = Quantity(
        type=str,
        shape=[],
        description='''
        Phonon DFPT solver method
        ''',
        a_legacy=LegacyDefinition(name='x_castep_DFPT_solver_method'))

    x_castep_phonon_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Phonon calculation tolerance (eV/A**2)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_phonon_tolerance'))

    x_castep_phonon_cycles = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Phonon calculation cycles
        ''',
        a_legacy=LegacyDefinition(name='x_castep_phonon_cycles'))

    x_castep_band_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Phonon band convergence tolerance window
        ''',
        a_legacy=LegacyDefinition(name='x_castep_band_tolerance'))


class x_castep_section_density_mixing_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_density_mixing_parameters'))

    x_castep_density_mixing_scheme = Quantity(
        type=str,
        shape=[],
        description='''
        density_mixing_scheme
        ''',
        a_legacy=LegacyDefinition(name='x_castep_density_mixing_scheme'))

    x_castep_density_mixing_length = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        density_mixing_scheme_length
        ''',
        a_legacy=LegacyDefinition(name='x_castep_density_mixing_length'))

    x_castep_charge_density_mixing_amplitude = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        charge_density_mixing_amplitude
        ''',
        a_legacy=LegacyDefinition(name='x_castep_charge_density_mixing_amplitude'))

    x_castep_cut_off_energy_for_mixing = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        charge_density_mixing_cut_off_energy_for_mixing (A)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cut_off_energy_for_mixing'))


class x_castep_section_population_analysis_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_population_analysis_parameters'))

    x_castep_population_analysis_cutoff = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Population_analysis_cutoff_(A)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_population_analysis_cutoff'))


class x_castep_section_ts_scf(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_scf'))

    x_castep_scf_ts_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        SCF_ts_energy Total
        ''',
        a_legacy=LegacyDefinition(name='x_castep_scf_ts_total'))

    x_castep_scf_ts_total_energy_free = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        SCF_ts_energy Total free
        ''',
        a_legacy=LegacyDefinition(name='x_castep_scf_ts_total_energy_free'))

    x_castep_scf_ts_T0 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        SCF_ts_energy T0 free
        ''',
        a_legacy=LegacyDefinition(name='x_castep_scf_ts_T0'))

    x_castep_section_ts_scf_iteration = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_scf_iteration'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_scf_iteration'))


class x_castep_section_k_band(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_k_band'))


class x_castep_section_k_points_1(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_k_points_1'))

    x_castep_store_k_points_1 = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing k points coordinates (fractional)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_k_points_1'))


class x_castep_section_k_points(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_k_points'))

    x_castep_store_k_points = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing k points coordinates (fractional)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_k_points'))


class x_castep_section_relativity_treatment(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_relativity_treatment'))

    x_castep_relativity_treatment_scf = Quantity(
        type=str,
        shape=[],
        description='''
        Relativity treatment in CASTEP convention
        ''',
        a_legacy=LegacyDefinition(name='x_castep_relativity_treatment_scf'))


class x_castep_section_scf_eigenvalues(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_scf_eigenvalues'))

    x_castep_store_scf_eigenvalues = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_scf_eigenvalues'))


class x_castep_section_SCF_iteration_frame(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_SCF_iteration_frame'))

    x_castep_frame_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        CASTEP_store_t_md_frame
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_time'))

    x_castep_SCF_frame_energy = Quantity(
        type=str,
        shape=[],
        description='''
        energy_frame_iterations
        ''',
        a_legacy=LegacyDefinition(name='x_castep_SCF_frame_energy'))

    x_castep_SCF_frame_energy_gain = Quantity(
        type=str,
        shape=[],
        description='''
        energy_frame_iterations_gain
        ''',
        a_legacy=LegacyDefinition(name='x_castep_SCF_frame_energy_gain'))

    x_castep_frame_time_scf_iteration_wall_end = Quantity(
        type=str,
        shape=[],
        description='''
        energy_frame_wall_end_time
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_time_scf_iteration_wall_end'))

    x_castep_frame_energy_free = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        energy_free
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_energy_free'))

    x_castep_frame_energy_total_T0 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        energy_free_corrected_for_finite_basis
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_energy_total_T0'))


class x_castep_section_scf_k_points(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_scf_k_points'))

    x_castep_store_scf_k_points = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_scf_k_points'))

    x_castep_section_scf_eigenvalues = SubSection(
        sub_section=SectionProxy('x_castep_section_scf_eigenvalues'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_scf_eigenvalues'))


class x_castep_section_spin_number(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_spin_number'))

    x_castep_spin_number = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_castep_spin_number'))


class x_castep_section_stress_tensor(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_stress_tensor'))

    x_castep_store_stress_tensor = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing stress tensor components
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_stress_tensor'))


class x_castep_section_time(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_time'))

    x_castep_calculation_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        castep_calculation_time
        ''',
        a_legacy=LegacyDefinition(name='x_castep_calculation_time'))

    x_castep_finalisation_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        castep_finalisation_time
        ''',
        a_legacy=LegacyDefinition(name='x_castep_finalisation_time'))

    x_castep_initialisation_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        x_castep_Initialisation_time
        ''',
        a_legacy=LegacyDefinition(name='x_castep_initialisation_time'))


class x_castep_section_raman_tensor(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_raman_tensor'))

    x_castep_store_raman_tensor = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing converged Raman susceptibility tensor
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_raman_tensor'))

    x_castep_raman_tensor = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        # TODO verify unit, original parser was ampere / amu
        unit='angstrom / unified_atomic_mass_unit',
        description='''
        Ramen tensor
        ''',
        a_legacy=LegacyDefinition(name='x_castep_raman_tensor'))


class x_castep_section_scf_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_scf_parameters'))

    x_castep_energy_threshold = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Energy Threshold
        ''',
        a_legacy=LegacyDefinition(name='x_castep_energy_threshold'))

    x_castep_max_iter = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of maximum iterations steps
        ''',
        a_legacy=LegacyDefinition(name='x_castep_max_iter'))

    x_castep_smearing_kind = Quantity(
        type=str,
        shape=[],
        description='''
        Smearing kind
        ''',
        a_legacy=LegacyDefinition(name='x_castep_smearing_kind'))

    x_castep_smearing_width = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Smearing width
        ''',
        a_legacy=LegacyDefinition(name='x_castep_smearing_width'))


class x_castep_section_electronic_spectroscpy_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_electronic_spectroscpy_parameters'))

    x_castep_theory_level = Quantity(
        type=str,
        shape=[],
        description='''
        Electronic spectroscopy parameters theory level
        ''',
        a_legacy=LegacyDefinition(name='x_castep_theory_level'))

    x_castep_spectroscopy_calculation = Quantity(
        type=str,
        shape=[],
        description='''
        Electronic spectroscopy parameters calculation type
        ''',
        a_legacy=LegacyDefinition(name='x_castep_spectroscopy_calculation'))

    x_castep_spec_max_iter = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Max number of iterations
        ''',
        a_legacy=LegacyDefinition(name='x_castep_spec_max_iter'))

    x_castep_spec_max_steps = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Max number of steps
        ''',
        a_legacy=LegacyDefinition(name='x_castep_spec_max_steps'))

    x_castep_spec_max_bands = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Max number of bands
        ''',
        a_legacy=LegacyDefinition(name='x_castep_spec_max_bands'))

    x_castep_spec_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Tolerance (eV)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_spec_tolerance'))


class x_castep_section_md(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_md'))

    x_castep_md_energies = Quantity(
        type=str,
        shape=[],
        description='''
        md_energy_components
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_energies'))

    x_castep_md_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        md_temp(K)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_temperature'))

    x_castep_md_pressure = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        md_pressure
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_pressure'))

    x_castep_md_cell_vectors = Quantity(
        type=str,
        shape=[],
        description='''
        md_cell_vectors
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_cell_vectors'))

    x_castep_md_cell_vectors_vel = Quantity(
        type=str,
        shape=[],
        description='''
        md_cell_vectors_velocities
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_cell_vectors_vel'))

    x_castep_md_stress_tensor = Quantity(
        type=str,
        shape=[],
        description='''
        md_stress_tensor
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_stress_tensor'))

    x_castep_md_positions = Quantity(
        type=str,
        shape=[],
        description='''
        md_positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_positions'))

    x_castep_md_forces = Quantity(
        type=str,
        shape=[],
        description='''
        md_forces
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_forces'))

    x_castep_md_lab = Quantity(
        type=str,
        shape=[],
        description='''
        md_lables
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_lab'))

    x_castep_md_veloc = Quantity(
        type=str,
        shape=[],
        description='''
        md_veloc
        ''',
        a_legacy=LegacyDefinition(name='x_castep_md_veloc'))


class x_castep_section_ts(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts'))

    x_castep_ts_energy_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_energy_components_total
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_energy_total'))

    x_castep_ts_cell_vectors = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        ts_cell_vectors
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_cell_vectors'))

    x_castep_ts_positions = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        ts_positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_positions'))

    x_castep_ts_forces = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        ts_forces
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_forces'))

    x_castep_ts_lab = Quantity(
        type=str,
        shape=[],
        description='''
        ts_lables
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_lab'))


class x_castep_section_ts_store(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_store'))

    x_castep_ts_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_energy_components_ts_store
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_energy'))

    x_castep_ts_cell_vectors_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_cell_vectors
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_cell_vectors_store'))

    x_castep_ts_positions_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_positions_store
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_positions_store'))

    x_castep_ts_forces_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_forces_store
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_forces_store'))


class x_castep_section_ts_final_store(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_final_store'))

    x_castep_ts_energy_final_store = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_energy_components_ts
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_energy_final_store'))

    x_castep_ts_forces_final_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_forces_final
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_forces_final_store'))

    x_castep_ts_positions_final_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_energy_positions_final
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_positions_final_store'))

    x_castep_ts_cell_vectors_final_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_cell_vectors_final
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_cell_vectors_final_store'))


class x_castep_section_ts_final(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_final'))

    x_castep_ts_energy_final = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_energy_final
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_energy_final'))

    x_castep_ts_cell_vectors_final = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        ts_cell_vectors
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_cell_vectors_final'))

    x_castep_ts_positions_final = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        ts_positions_final
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_positions_final'))

    x_castep_ts_forces_final = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        ts_force_finals
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_forces_final'))


class x_castep_section_ts_product(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_product'))

    x_castep_ts_energy_product = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_energy_prod
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_energy_product'))

    x_castep_ts_forces_product = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        ts_energy_forces_pro
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_forces_product'))

    x_castep_ts_positions_product = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        ts_energy_positions_pro
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_positions_product'))

    x_castep_ts_cell_vectors_product = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        ts_product_cell_vectors
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_cell_vectors_product'))


class x_castep_section_ts_product_store(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_product_store'))

    x_castep_ts_energy_product_store = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_energy_components_pro
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_energy_product_store'))

    x_castep_ts_forces_pro_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_energy_forces_pro_store
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_forces_pro_store'))

    x_castep_ts_positions_pro_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_energy_positions_pro_store
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_positions_pro_store'))

    x_castep_ts_cell_vectors_pro_store = Quantity(
        type=str,
        shape=[],
        description='''
        ts_cell_vectors_pro
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_cell_vectors_pro_store'))


class x_castep_section_ts_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_ts_parameters'))

    x_castep_ts_method = Quantity(
        type=str,
        shape=[],
        description='''
        ts_method
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_method'))

    x_castep_ts_protocol = Quantity(
        type=str,
        shape=[],
        description='''
        ts_protocol
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_protocol'))

    x_castep_ts_number_qst = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_qst_iterations
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_number_qst'))

    x_castep_ts_number_cg = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_number_of_cg_iterations
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_number_cg'))

    x_castep_ts_force_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_force_tolerance (eV/A)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_force_tolerance'))

    x_castep_ts_displacement_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_displacement_tolerance (A)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_displacement_tolerance'))


class x_castep_section_DFT_SEDC(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_DFT_SEDC'))

    x_castep_correction_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        correlation energy
        ''',
        a_legacy=LegacyDefinition(name='x_castep_correction_energy'))

    x_castep_de_atom = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        dE/atom
        ''',
        a_legacy=LegacyDefinition(name='x_castep_de_atom'))

    x_castep_dfmax_atom = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        dfmax/atom
        ''',
        a_legacy=LegacyDefinition(name='x_castep_dfmax_atom'))

    x_castep_structure_energy_corr = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        structure energy correction
        ''',
        a_legacy=LegacyDefinition(name='x_castep_structure_energy_corr'))

    x_castep_PBC_image_inter_corr = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        PBC image interaction corr.
        ''',
        a_legacy=LegacyDefinition(name='x_castep_PBC_image_inter_corr'))

    x_castep_total_energy_correction = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        total energy correction
        ''',
        a_legacy=LegacyDefinition(name='x_castep_total_energy_correction'))

    x_castep_total_fmax_correction = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        correction F max ev/A
        ''',
        a_legacy=LegacyDefinition(name='x_castep_total_fmax_correction'))

    x_castep_shell = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        shell
        ''',
        a_legacy=LegacyDefinition(name='x_castep_shell'))

    x_castep_total_dispersion_corrected_free_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        total_dispersion_corrected_free_energy
        ''',
        a_legacy=LegacyDefinition(name='x_castep_total_dispersion_corrected_free_energy'))

    x_castep_disp_corrected_energy_total_T0 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        dispersion corrected zero point
        ''',
        a_legacy=LegacyDefinition(name='x_castep_disp_corrected_energy_total_T0'))


class x_castep_section_van_der_Waals_parameters(MSection):
    '''
    -
    '''

    m_def = Section(validate=False, a_legacy=LegacyDefinition(name='x_castep_section_van_der_Waals_parameters'))

    x_castep_disp_method_name = Quantity(
        type=str,
        shape=[],
        description='''
        Name type
        ''',
        a_legacy=LegacyDefinition(name='x_castep_disp_method_name'))

    x_castep_Parameter_d = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Parameter for dispersion method G06
        ''',
        a_legacy=LegacyDefinition(name='x_castep_Parameter_d'))

    x_castep_Parameter_LAMBDA = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Parameter for dispersion method OBS
        ''',
        a_legacy=LegacyDefinition(name='x_castep_Parameter_LAMBDA'))

    x_castep_Parameter_n = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Parameter for dispersion method OBS
        ''',
        a_legacy=LegacyDefinition(name='x_castep_Parameter_n'))

    x_castep_Parameter_s6 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Parameter for dispersion method G06
        ''',
        a_legacy=LegacyDefinition(name='x_castep_Parameter_s6'))

    x_castep_Parameter_sR = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Parameter for dispersion method TS
        ''',
        a_legacy=LegacyDefinition(name='x_castep_Parameter_sR'))


class section_run(public.section_run):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_run'))

    x_castep_atom_forces = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        Forces on the atoms as minus gradient of energy_total, including forces' unitary-
        transformation (rigid body) filtering and including constraints, if present. The
        derivatives with respect to displacements of the nuclei in the gradient are
        evaluated according to the coordinate system defined in coordinate_system. In
        addition, these forces are obtained by filtering out the unitary transformations
        (translations of the center of mass and rigid rotations of the whole system, when
        non periodic), atom_forces_raw for the unfiltered counterpart. Furthermore, forces
        due to constraints like fixed atoms, distances, angles, dihedrals, and so on, are
        here included (see atom_forces_raw for the unfiltered counterpart).
        ''',
        a_legacy=LegacyDefinition(name='x_castep_atom_forces'))

    x_castep_basis_set_planewave_cutoff_iteration_0 = Quantity(
        type=str,
        shape=[],
        description='''
        cutoff at iteration 0 of geometry optimisation
        ''',
        a_legacy=LegacyDefinition(name='x_castep_basis_set_planewave_cutoff_iteration_0'))

    x_castep_compiler = Quantity(
        type=str,
        shape=[],
        description='''
        Compiler name
        ''',
        a_legacy=LegacyDefinition(name='x_castep_compiler'))

    x_castep_constants_reference = Quantity(
        type=str,
        shape=[],
        description='''
        Fundamental constant data source
        ''',
        a_legacy=LegacyDefinition(name='x_castep_constants_reference'))

    x_castep_fft_library = Quantity(
        type=str,
        shape=[],
        description='''
        fft library name
        ''',
        a_legacy=LegacyDefinition(name='x_castep_fft_library'))

    x_castep_geom_converged = Quantity(
        type=str,
        shape=[],
        description='''
        CASTEP_geom_converged
        ''',
        a_legacy=LegacyDefinition(name='x_castep_geom_converged'))

    x_castep_maths_library = Quantity(
        type=str,
        shape=[],
        description='''
        Maths library name
        ''',
        a_legacy=LegacyDefinition(name='x_castep_maths_library'))

    x_castep_program_compilation_date = Quantity(
        type=str,
        shape=[],
        description='''
        Compilation date (string)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_program_compilation_date'))

    x_castep_program_compilation_time = Quantity(
        type=str,
        shape=[],
        description='''
        Compilation time (string)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_program_compilation_time'))

    x_castep_program_execution_date = Quantity(
        type=str,
        shape=[],
        description='''
        Run execution date (string)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_program_execution_date'))

    x_castep_program_execution_time = Quantity(
        type=str,
        shape=[],
        description='''
        Run execution start time (string)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_program_execution_time'))

    x_castep_crystal_point_group = Quantity(
        type=str,
        shape=[],
        description='''
        Point group of the crystal (Schoenflies notation)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_crystal_point_group'))

    x_castep_space_group = Quantity(
        type=str,
        shape=[],
        description='''
        Point space of the crystal
        ''',
        a_legacy=LegacyDefinition(name='x_castep_space_group'))

    x_castep_ts_path = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_path_number
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_path'))

    x_castep_ts_path_ts_final = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_path_final
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_path_ts_final'))

    x_castep_ts_path_product = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        ts_path_pro
        ''',
        a_legacy=LegacyDefinition(name='x_castep_ts_path_product'))

    x_castep_store_atom_forces_band = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing converged atom forces (ev/A)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_forces_band'))

    x_castep_elec_methd = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing electronic structure method
        ''',
        a_legacy=LegacyDefinition(name='x_castep_elec_methd'))

    x_castep_section_vibrational_frequencies = SubSection(
        sub_section=SectionProxy('x_castep_section_vibrational_frequencies'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_vibrational_frequencies'))

    x_castep_section_band_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_band_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_band_parameters'))

    x_castep_section_core_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_core_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_core_parameters'))

    x_castep_section_collect_scf_eigenvalues = SubSection(
        sub_section=SectionProxy('x_castep_section_collect_scf_eigenvalues'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_collect_scf_eigenvalues'))

    x_castep_section_population_analysis = SubSection(
        sub_section=SectionProxy('x_castep_section_population_analysis'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_population_analysis'))

    x_castep_section_geom_optimisation_method = SubSection(
        sub_section=SectionProxy('x_castep_section_geom_optimisation_method'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_geom_optimisation_method'))

    x_castep_section_optics_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_optics_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_optics_parameters'))

    x_castep_section_tddft_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_tddft_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_tddft_parameters'))

    x_castep_section_phonons = SubSection(
        sub_section=SectionProxy('x_castep_section_phonons'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_phonons'))

    x_castep_section_density_mixing_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_density_mixing_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_density_mixing_parameters'))

    x_castep_section_population_analysis_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_population_analysis_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_population_analysis_parameters'))

    x_castep_section_SCF_iteration_frame = SubSection(
        sub_section=SectionProxy('x_castep_section_SCF_iteration_frame'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_SCF_iteration_frame'))

    x_castep_section_time = SubSection(
        sub_section=SectionProxy('x_castep_section_time'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_time'))

    x_castep_section_raman_tensor = SubSection(
        sub_section=SectionProxy('x_castep_section_raman_tensor'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_raman_tensor'))

    x_castep_section_scf_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_scf_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_scf_parameters'))

    x_castep_section_electronic_spectroscpy_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_electronic_spectroscpy_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_electronic_spectroscpy_parameters'))

    x_castep_section_md = SubSection(
        sub_section=SectionProxy('x_castep_section_md'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_md'))

    x_castep_section_ts = SubSection(
        sub_section=SectionProxy('x_castep_section_ts'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts'))

    x_castep_section_ts_store = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_store'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_store'))

    x_castep_section_ts_final_store = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_final_store'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_final_store'))

    x_castep_section_ts_final = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_final'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_final'))

    x_castep_section_ts_product = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_product'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_product'))

    x_castep_section_ts_product_store = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_product_store'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_product_store'))

    x_castep_section_ts_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_parameters'))

    x_castep_section_van_der_Waals_parameters = SubSection(
        sub_section=SectionProxy('x_castep_section_van_der_Waals_parameters'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_van_der_Waals_parameters'))


class section_system(public.section_system):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_system'))

    x_castep_net_charge = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Net charge of system
        ''',
        a_legacy=LegacyDefinition(name='x_castep_net_charge'))

    x_castep_number_of_bands = Quantity(
        type=np.dtype(np.int32),
        shape=[],
        description='''
        Number of bands
        ''',
        a_legacy=LegacyDefinition(name='x_castep_number_of_bands'))

    x_castep_number_of_electrons = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number of electrons
        ''',
        a_legacy=LegacyDefinition(name='x_castep_number_of_electrons'))

    x_castep_atom_positions = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        Storing atomic positions in fractional coordinates
        ''',
        a_legacy=LegacyDefinition(name='x_castep_atom_positions'))

    x_castep_cell_volume = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        CASTEP_cell_volume
        ''',
        a_legacy=LegacyDefinition(name='x_castep_cell_volume'))

    x_castep_optimised_atom_labels = Quantity(
        type=str,
        shape=['number_of_atoms'],
        description='''
        Temporary storing atomic positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_optimised_atom_labels'))

    x_castep_optimised_atom_positions = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        Storing atomic optimised positions in fractional coordinates
        ''',
        a_legacy=LegacyDefinition(name='x_castep_optimised_atom_positions'))

    x_castep_velocities_cell_vector = Quantity(
        type=np.dtype(np.float64),
        shape=[3, 3],
        description='''
        cell vector velocities
        ''',
        a_legacy=LegacyDefinition(name='x_castep_velocities_cell_vector'))

    x_castep_store_atom_labels = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atom labels
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_labels'))

    x_castep_store_atom_number = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atom labels
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_number'))

    x_castep_store_atom_positions = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atomic positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_positions'))

    x_castep_store_atom_ionic_velocities = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atomic positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_ionic_velocities'))

    x_castep_atom_ionic_velocities = Quantity(
        type=np.dtype(np.float64),
        shape=['number_of_atoms', 3],
        description='''
        Temporary storing atomic positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_atom_ionic_velocities'))

    x_castep_store_optimised_atom_labels = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atomic positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_optimised_atom_labels'))

    x_castep_store_optimised_atom_positions = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atomic positions
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_optimised_atom_positions'))

    x_castep_section_tddft = SubSection(
        sub_section=SectionProxy('x_castep_section_tddft'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_tddft'))

    x_castep_section_atom_ionic_velocities = SubSection(
        sub_section=SectionProxy('x_castep_section_atom_ionic_velocities'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_atom_ionic_velocities'))

    x_castep_section_atom_positions_optim = SubSection(
        sub_section=SectionProxy('x_castep_section_atom_positions_optim'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_atom_positions_optim'))

    x_castep_section_atom_positions = SubSection(
        sub_section=SectionProxy('x_castep_section_atom_positions'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_atom_positions'))

    x_castep_section_cell_optim = SubSection(
        sub_section=SectionProxy('x_castep_section_cell_optim'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_cell_optim'))

    x_castep_section_cell = SubSection(
        sub_section=SectionProxy('x_castep_section_cell'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_cell'))


class section_basis_set_cell_dependent(public.section_basis_set_cell_dependent):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_basis_set_cell_dependent'))

    x_castep_basis_set_planewave_cutoff = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing plane wave cutoff as string
        ''',
        a_legacy=LegacyDefinition(name='x_castep_basis_set_planewave_cutoff'))

    x_castep_size_std_grid = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        size of standard grid (eV)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_size_std_grid'))

    x_castep_size_fine_grid = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        sise of fine grid (1/A)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_size_fine_grid'))


class section_single_configuration_calculation(public.section_single_configuration_calculation):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_single_configuration_calculation'))

    x_castep_enthalpy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        CASTEP_Enthalpy
        ''',
        a_legacy=LegacyDefinition(name='x_castep_enthalpy'))

    x_castep_frequency = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        CASTEP_frequency (cm-1)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frequency'))

    x_castep_improved_energy_total = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        md_forces
        ''',
        a_legacy=LegacyDefinition(name='x_castep_improved_energy_total'))

    x_castep_frame_time_0 = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        CASTEP_store_t_md_frame
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_time_0'))

    x_castep_geom_iteration_index = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Index for number of iterations in geometry optimisation
        ''',
        a_legacy=LegacyDefinition(name='x_castep_geom_iteration_index'))

    x_castep_store_atom_forces = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing converged atom forces (ev/A)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_forces'))

    x_castep_initial_scf_iteration_wall_time = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Initial SCF iteration wall time
        ''',
        a_legacy=LegacyDefinition(name='x_castep_initial_scf_iteration_wall_time'))

    x_castep_total_dispersion_corrected_energy = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Total electronic energy that includes dispersion energy computed with
        Disp_method_name not corrected for finite basis-set
        ''',
        a_legacy=LegacyDefinition(name='x_castep_total_dispersion_corrected_energy'))

    x_castep_total_energy_corrected_for_finite_basis = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        CASTEP_total_energy_corrected_for_finite_basis
        ''',
        a_legacy=LegacyDefinition(name='x_castep_total_energy_corrected_for_finite_basis'))

    x_castep_total_energy_corrected_for_finite_basis_store = Quantity(
        type=str,
        shape=[],
        description='''
        CASTEP_total_energy_corrected_for_finite_basis_store
        ''',
        a_legacy=LegacyDefinition(name='x_castep_total_energy_corrected_for_finite_basis_store'))

    energy_total_scf_iteration_list = Quantity(
        type=np.dtype(np.float64),
        shape=[-1],
        description='''
        Total electronic energy calculated with XC_method_scf during the scf iterations is
        stored in a list
        ''',
        a_legacy=LegacyDefinition(name='energy_total_scf_iteration_list'))

    x_castep_section_ts_scf = SubSection(
        sub_section=SectionProxy('x_castep_section_ts_scf'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_ts_scf'))

    x_castep_section_scf_k_points = SubSection(
        sub_section=SectionProxy('x_castep_section_scf_k_points'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_scf_k_points'))

    x_castep_section_spin_number = SubSection(
        sub_section=SectionProxy('x_castep_section_spin_number'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_spin_number'))

    x_castep_section_stress_tensor = SubSection(
        sub_section=SectionProxy('x_castep_section_stress_tensor'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_stress_tensor'))

    x_castep_section_DFT_SEDC = SubSection(
        sub_section=SectionProxy('x_castep_section_DFT_SEDC'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_DFT_SEDC'))


class section_method(public.section_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_method'))

    x_castep_functional_and_weight = Quantity(
        type=str,
        shape=[],
        description='''
        XC functional+weight in CASTEP convention
        ''',
        a_legacy=LegacyDefinition(name='x_castep_functional_and_weight'))

    x_castep_section_functional_definition = SubSection(
        sub_section=SectionProxy('x_castep_section_functional_definition'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_functional_definition'))

    x_castep_section_functionals = SubSection(
        sub_section=SectionProxy('x_castep_section_functionals'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_functionals'))

    x_castep_section_relativity_treatment = SubSection(
        sub_section=SectionProxy('x_castep_section_relativity_treatment'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_relativity_treatment'))


class section_sampling_method(public.section_sampling_method):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_sampling_method'))

    x_castep_geometry_stress_com_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='pascal',
        description='''
        tolerance for stress components in geometry optimisation (GPa)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_geometry_stress_com_tolerance'))

    x_castep_thermostat_target_temperature = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='kelvin',
        description='''
        thermostat_target_temperature(K)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_thermostat_target_temperature'))

    x_castep_barostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        barostat_type
        ''',
        a_legacy=LegacyDefinition(name='x_castep_barostat_type'))

    x_castep_thermostat_type = Quantity(
        type=str,
        shape=[],
        description='''
        thermostat_type
        ''',
        a_legacy=LegacyDefinition(name='x_castep_thermostat_type'))

    x_castep_thermostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        thermostat_type
        ''',
        a_legacy=LegacyDefinition(name='x_castep_thermostat_tau'))

    x_castep_barostat_tau = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        barostat_tau
        ''',
        a_legacy=LegacyDefinition(name='x_castep_barostat_tau'))

    x_castep_integrator_dt = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='second',
        description='''
        MD_time_step (ps)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_integrator_dt'))

    x_castep_number_of_steps_requested = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        MD_time_step_number
        ''',
        a_legacy=LegacyDefinition(name='x_castep_number_of_steps_requested'))

    x_castep_max_number_of_steps = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        Number_of iterations geom_optim
        ''',
        a_legacy=LegacyDefinition(name='x_castep_max_number_of_steps'))

    x_castep_frame_pressure = Quantity(
        type=str,
        shape=[],
        description='''
        MD_pressure
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_pressure'))

    x_castep_frame_energy_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        MD_scf_energy tolerance (eV)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_energy_tolerance'))

    x_castep_frame_eigen_tolerance = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        unit='joule',
        description='''
        MD_scf_eigen tolerance (eV)
        ''',
        a_legacy=LegacyDefinition(name='x_castep_frame_eigen_tolerance'))


class section_k_band(public.section_k_band):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_k_band'))

    x_castep_k_path = Quantity(
        type=np.dtype(np.float64),
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_castep_k_path'))

    x_castep_store_k_path = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_k_path'))

    x_castep_store_k_label = Quantity(
        type=str,
        shape=[],
        description='''
        -
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_k_label'))

    x_castep_section_eigenvalues_1 = SubSection(
        sub_section=SectionProxy('x_castep_section_eigenvalues_1'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_eigenvalues_1'))

    x_castep_section_eigenvalues = SubSection(
        sub_section=SectionProxy('x_castep_section_eigenvalues'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_eigenvalues'))

    x_castep_section_k_band = SubSection(
        sub_section=SectionProxy('x_castep_section_k_band'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_k_band'))

    x_castep_section_k_points_1 = SubSection(
        sub_section=SectionProxy('x_castep_section_k_points_1'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_k_points_1'))

    x_castep_section_k_points = SubSection(
        sub_section=SectionProxy('x_castep_section_k_points'),
        repeats=True,
        a_legacy=LegacyDefinition(name='x_castep_section_k_points'))


class section_topology(common.section_topology):

    m_def = Section(validate=False, extends_base_section=True, a_legacy=LegacyDefinition(name='section_topology'))

    x_castep_store_atom_mass = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atom mass
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_mass'))

    x_castep_store_atom_name = Quantity(
        type=str,
        shape=[],
        description='''
        Temporary storing atom name
        ''',
        a_legacy=LegacyDefinition(name='x_castep_store_atom_name'))


m_package.__init_metainfo__()
