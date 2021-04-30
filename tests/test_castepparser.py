#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD. See https://nomad-lab.eu for further info.
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

import pytest
import numpy as np

from nomad.datamodel import EntryArchive
from castepparser.castep_parser import CastepParser


def approx(value, abs=0, rel=1e-6):
    return pytest.approx(value, abs=abs, rel=rel)


@pytest.fixture(scope='module')
def parser():
    return CastepParser()


def test_single_point(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Si8.castep', archive, None)

    sec_run = archive.section_run[0]
    assert sec_run.program_version == '16.1'
    assert sec_run.x_castep_constants_reference == 'CODATA 2010'
    assert sec_run.time_run_date_start.magnitude == 1455286325.0
    assert sec_run.section_basis_set_cell_dependent[0].basis_set_cell_dependent_name == 'PW_7'

    sec_method = sec_run.section_method[0]
    assert sec_method.number_of_spin_channels == 1
    assert sec_method.section_XC_functionals[1].XC_functional_name == 'GGA_C_PBE'
    assert sec_method.smearing_kind == 'gaussian'

    assert sec_run.section_sampling_method[0].sampling_method == 'single_point'

    sec_scc = sec_run.section_single_configuration_calculation[0]
    assert np.shape(sec_scc.atom_forces) == (8, 3)
    assert np.count_nonzero(sec_scc.atom_forces) == 0
    assert sec_scc.energy_total.magnitude == approx(-2.21372403e-16)
    sec_scfs = sec_scc.section_scf_iteration
    assert len(sec_scfs) == 12
    assert sec_scfs[3].energy_total_scf_iteration.magnitude == approx(-2.21530505e-16)
    assert sec_scfs[8].energy_reference_fermi_iteration[0].magnitude == approx(8.53007633e-19)
    assert sec_scfs[6].time_scf_iteration.magnitude == 12.70

    sec_system = sec_run.section_system[0]
    assert sec_system.atom_positions[2][1].magnitude == approx(2.715e-10)
    assert sec_system.x_castep_section_atom_positions[0].x_castep_cell_length_a == approx(5.43e-10)

    assert len(sec_run.section_topology[0].section_atom_type) == 1

    sec_mulliken = sec_run.x_castep_section_population_analysis
    assert len(sec_mulliken) == 8
    assert sec_mulliken[4].x_castep_orbital_p == 2.66

    assert sec_run.x_castep_section_density_mixing_parameters[0].x_castep_density_mixing_length == 20
    assert sec_run.x_castep_section_time[0].x_castep_finalisation_time == 0.01


def test_dmd(parser):
    archive = EntryArchive()
    parser.parse('tests/data/TiO2-geom.castep', archive, None)

    sec_sampling = archive.section_run[0].section_sampling_method[0]
    assert sec_sampling.sampling_method == 'geometry_optimization'
    assert sec_sampling.geometry_optimization_method == 'damped MD'
    assert sec_sampling.geometry_optimization_threshold_force.magnitude == approx(8.01088317e-11)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 23
    assert sec_sccs[7].energy_total.magnitude == approx(-3.67496146e-16)
    assert len(sec_sccs[20].section_scf_iteration) == 13
    assert sec_sccs[17].section_scf_iteration[3].energy_total_scf_iteration.magnitude == approx(-3.67497215e-16)
    assert sec_sccs[3].atom_forces[0][1].magnitude == approx(-4.49314415e-10,)
    assert sec_sccs[22].energy_total.magnitude == approx(-3.67497218e-16)

    sec_systems = archive.section_run[0].section_system
    assert sec_systems[12].atom_positions[3][0].magnitude == approx(9.074282e-11)
    assert sec_systems[0].x_castep_number_of_electrons == 32


def test_md(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Si8-md-NPT.castep', archive, None)

    sec_sampling = archive.section_run[0].section_sampling_method[0]
    assert sec_sampling.sampling_method == 'molecular_dynamics'
    assert sec_sampling.ensemble_type == 'NPT'
    assert sec_sampling.x_castep_thermostat_type == 'Nose-Hoover chain thermostat'
    assert sec_sampling.x_castep_frame_energy_tolerance.magnitude == approx(1.60217663e-24)

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 13
    assert sec_sccs[2].energy_total.magnitude == approx(-1.37059981e-16)
    assert sec_sccs[6].stress_tensor[2][1].magnitude == approx(4.06626e+09)
    assert sec_sccs[9].pressure.magnitude == approx(1.111e+09)
    assert sec_sccs[11].energy_total_T0.magnitude == approx(-1.37069057e-16)
    assert len(sec_sccs[12].section_scf_iteration) == 7
    assert sec_sccs[7].section_scf_iteration[3].energy_change_scf_iteration.magnitude == approx(-1.90981043e-21)

    sec_systems = archive.section_run[0].section_system
    assert len(sec_systems) == 13
    assert sec_systems[0].atom_positions[7][2].magnitude == approx(1.3575e-10)
    assert sec_systems[1].atom_velocities[1][1].magnitude == approx(324.536)
    assert sec_systems[12].lattice_vectors[2][2].magnitude == approx(5.5475876e-10)


def test_eigenvalues(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Fe.castep', archive, None)

    sec_eigenvalues = archive.section_run[0].section_single_configuration_calculation[0].section_eigenvalues[0]
    assert np.shape(sec_eigenvalues.eigenvalues_values) == (2, 118, 6)
    assert sec_eigenvalues.eigenvalues_values[1][38][4].magnitude == approx(1.30819997e-18)
    assert sec_eigenvalues.eigenvalues_kpoints[22][1] == 0.289474


def test_bandstructure(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Dispersions/Si2.castep', archive, None)

    sec_band_segment = archive.section_run[0].section_single_configuration_calculation[0].section_k_band[0].section_k_band_segment
    assert len(sec_band_segment) == 5
    assert sec_band_segment[3].band_segm_labels == ['X', 'W']
    assert sec_band_segment[1].band_energies[0][-1][12].magnitude == approx(2.17418526e-18)
    assert sec_band_segment[4].band_k_points[2][1] == 0.300000


def test_vibration(parser):
    archive = EntryArchive()
    parser.parse('tests/data/BC2N-Pmm2-Raman.castep', archive, None)

    sec_vibration = archive.section_run[0].x_castep_section_vibrational_frequencies
    assert len(sec_vibration) == 2
    assert sec_vibration[1].x_castep_vibrational_frequencies[2] == approx(0.461821)
    assert sec_vibration[0].x_castep_raman_activity[3] == approx(21.0162567)
    assert sec_vibration[1].x_castep_ir_intensity[6] == approx(2.7078705)
    assert sec_vibration[0].x_castep_raman_active[10] == 'Y'

    sec_raman = archive.section_run[0].x_castep_section_raman_tensor
    assert len(sec_raman) == 12
    assert sec_raman[9].x_castep_raman_tensor[2][0].magnitude == approx(0.1834 * 0.5)


def test_tss(parser):
    archive = EntryArchive()
    parser.parse('tests/data/h2-lst.castep', archive, None)

    assert archive.section_run[0].section_sampling_method[0].sampling_method == 'geometry_optimization'

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 27
    assert sec_sccs[20].energy_total.magnitude == approx(-4.72809265e-18)


def test_bfgs(parser):
    archive = EntryArchive()
    parser.parse('tests/data/Si2_opt.castep', archive, None)

    assert archive.section_run[0].section_sampling_method[0].sampling_method == 'geometry_optimization'

    sec_sccs = archive.section_run[0].section_single_configuration_calculation
    assert len(sec_sccs) == 9
    sec_sccs[7].pressure.magnitude == approx(400000.0)


def test_di(parser):
    # TODO implement test cannot find an example
    pass
