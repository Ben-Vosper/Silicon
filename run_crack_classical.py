"""
run_crack_classical.py

Script to run classical molecular dynamics for a crack slab,
incrementing the load in small steps until fracture starts.

James Kermode <james.kermode@kcl.ac.uk>
January 2013
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
                         
# ******* Start of parameters ***********

input_file = 'crack.xyz'         # File from which to read crack slab structure
sim_T = 300.0*units.kB           # Simulation temperature
nsteps = 10000                   # Total number of timesteps to run for
timestep = 1.0*units.fs          # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang      # Amount by which potential cutoff is increased
                                 # for neighbour calculations
tip_move_tol = 10.0              # Distance tip has to move before crack 
                                 # is taken to be running
strain_rate = 1e-5*(1/units.fs)  # Strain rate
traj_file = 'traj.nc'            # Trajectory output file (NetCDF format)
traj_interval = 10               # Number of time steps between
                                 # writing output frames


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

mm_pot = Potential("IP SW", cutoff_skin=cutoff_skin)
mm_pot.set_default_properties(['stresses'])

atoms.set_calculator(mm_pot)

# ********* Setup and run MD ***********

MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)

dynamics = VelocityVerlet(atoms, timestep)

def printstatus():
    if dynamics.nsteps == 1:
        print """
State      Time/fs    Temp/K     Strain      G/(J/m^2)  CrackPos/A D(CrackPos)/A 
---------------------------------------------------------------------------------"""

    log_format = ('%(label)-4s%(time)12.1f%(temperature)12.6f'+
                  '%(strain)12.5f%(G)12.4f%(crack_pos_x)12.2f    (%(d_crack_pos_x)+5.2f)')

    atoms.info['label'] = 'D'                  # Label for the status line
    atoms.info['time'] = dynamics.get_time()/units.fs
    atoms.info['temperature'] = (atoms.get_kinetic_energy() /
                                 (1.5*units.kB*len(atoms)))
    atoms.info['strain'] = get_strain(atoms)
    atoms.info['G'] = get_energy_release_rate(atoms)/(units.J/units.m**2)
    
    crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)
    atoms.info['crack_pos_x'] = crack_pos[0]
    atoms.info['d_crack_pos_x'] = crack_pos[0] - orig_crack_pos[0]

    print log_format % atoms.info
dynamics.attach(printstatus)

def check_if_cracked(atoms):
    crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot)

    # stop straining if crack has advanced more than tip_move_tol
    if not atoms.info['is_cracked'] and (crack_pos[0] - orig_crack_pos[0]) > tip_move_tol:
        atoms.info['is_cracked'] = True
        del atoms.constraints[atoms.constraints.index(strain_atoms)]
dynamics.attach(check_if_cracked, 1, atoms)

trajectory = AtomsWriter(traj_file)
dynamics.attach(trajectory, traj_interval, atoms)

dynamics.run(nsteps)
