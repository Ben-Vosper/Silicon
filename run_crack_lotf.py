"""
run_crack_lotf_3.py

Script to run LOTF molecular dynamics for a crack slab,
incrementing the load in small steps until fracture starts.

This version checks the predictor-corrector force errors.

James Kermode <james.kermode@kcl.ac.uk>
February 2013
"""

import numpy as np

from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import ase.units as units

from quippy import set_fortran_indexing
from quippy.atoms import Atoms
from quippy.potential import Potential
from quippy.io import AtomsWriter

from quippy.crack import (get_strain,
                          get_energy_release_rate,
                          ConstantStrainRate,
                          find_crack_tip_stress_field)
from quippy.potential import ForceMixingPotential
from quippy.lotf import LOTFDynamics, update_hysteretic_qm_region

# ******* Start of parameters ***********

input_file = 'crack.xyz'         # File from which to read crack slab structure
sim_T = 300.0*units.kB           # Simulation temperature
nsteps = 200                     # Total number of timesteps to run for (reduced)
timestep = 1.0*units.fs          # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang      # Amount by which potential cutoff is increased
                                 # for neighbour calculations
tip_move_tol = 10.0              # Distance tip has to move before crack 
                                 # is taken to be running
strain_rate = 1e-5*(1/units.fs)  # Strain rate
traj_file = 'traj_q.nc'            # Trajectory output file (NetCDF format)
traj_interval = 10               # Number of time steps between
                                 # writing output frames

qm_init_args = 'TB DFTB'         # Initialisation arguments for QM potential
qm_inner_radius = 8.0*units.Ang  # Inner hysteretic radius for QM region
qm_outer_radius = 10.0*units.Ang # Inner hysteretic radius for QM region
extrapolate_steps = 10           # Number of steps for predictor-corrector
                                 # interpolation and extrapolation

set_fortran_indexing(False)

atoms = Atoms(input_file)

orig_height = atoms.info['OrigHeight']
orig_crack_pos = atoms.info['CrackPos'].copy()

top = atoms.positions[:, 1].max()
bottom = atoms.positions[:, 1].min()
left = atoms.positions[:, 0].min()
right = atoms.positions[:, 0].max()

fixed_mask = ((abs(atoms.positions[:, 1] - top) < 1.0) |
              (abs(atoms.positions[:, 1] - bottom) < 1.0))
fix_atoms = FixAtoms(mask=fixed_mask)

strain_atoms = ConstantStrainRate(orig_height, strain_rate*timestep)
atoms.set_constraint([fix_atoms, strain_atoms])

mm_pot = Potential("IP SW", param_filename="params.xml", cutoff_skin=cutoff_skin)
mm_pot.set_default_properties(['stresses'])

# Density functional tight binding (DFTB) potential
qm_pot = Potential(qm_init_args, param_filename="params.xml")
qm_list = update_hysteretic_qm_region(atoms, [], orig_crack_pos, qm_inner_radius, qm_outer_radius)
qmmm_pot = ForceMixingPotential(pot1=mm_pot,
                                pot2=qm_pot)
atoms.set_calculator(qmmm_pot)
qmmm_pot.set_qm_atoms(qm_list)


# *** Set up the initial QM region ****


#
#
# MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)
#
# dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps, check_force_error=True)
#
# def printstatus():
#     if dynamics.nsteps == 1:
#         print """
# State      Time/fs    Temp/K     Strain      G/(J/m^2)  CrackPos/A D(CrackPos)/A
# ---------------------------------------------------------------------------------"""
#
#     log_format = ('%(label)-4s%(time)12.1f%(temperature)12.6f'+
#                   '%(strain)12.5f%(G)12.4f%(crack_pos_x)12.2f    (%(d_crack_pos_x)+5.2f)')
#
#     atoms.info['label'] = dynamics.state_label  # Label for the status line
#     atoms.info['time'] = dynamics.get_time()/units.fs
#     atoms.info['temperature'] = (atoms.get_kinetic_energy() /
#                                  (1.5*units.kB*len(atoms)))
#     atoms.info['strain'] = get_strain(atoms)
#     atoms.info['G'] = get_energy_release_rate(atoms)/(units.J/units.m**2)
#
#     crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)
#     atoms.info['crack_pos_x'] = crack_pos[0]
#     atoms.info['d_crack_pos_x'] = crack_pos[0] - orig_crack_pos[0]
#
#     print log_format % atoms.info
# dynamics.attach(printstatus)
#
# def check_if_cracked(atoms):
#     crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)
#
#     # stop straining if crack has advanced more than tip_move_tol
#     if not atoms.info['is_cracked'] and (crack_pos[0] - orig_crack_pos[0]) > tip_move_tol:
#         atoms.info['is_cracked'] = True
#         del atoms.constraints[atoms.constraints.index(strain_atoms)]
# dynamics.attach(check_if_cracked, 1, atoms)
#
# def update_qm_region(atoms):
#    crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)
#    qm_list = qmmm_pot.get_qm_atoms()
#    qm_list = update_hysteretic_qm_region(atoms, qm_list, crack_pos,
#                                          qm_inner_radius, qm_outer_radius)
#    qmmm_pot.set_qm_atoms(qm_list)
# dynamics.set_qm_update_func(update_qm_region)
#
# trajectory = AtomsWriter(traj_file)
#
# def traj_writer(dynamics):
#    if dynamics.state == LOTFDynamics.Interpolation:
#       trajectory.write(dynamics.atoms)
# dynamics.attach(traj_writer, traj_interval, dynamics)
# #
# # def log_pred_corr_errors(dynamics, logfile):
# #     logfile.write('%s err %10.1f%12.6f%12.6f\n' % (dynamics.state_label,
# #                                                    dynamics.get_time()/units.fs,
# #                                                    dynamics.rms_force_error,
# #                                                    dynamics.max_force_error))
# # logfile = open('pred-corr-error.txt', 'w')
# # dynamics.attach(log_pred_corr_errors, 1, dynamics, logfile)
# #
# # # Start running!
# dynamics.run(nsteps)
#
# #logfile.close()
