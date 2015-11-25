"""
make_crack.py

Script to generate a crack slab, and apply initial strain ramp

James Kermode <james.kermode@kcl.ac.uk>
January 2013
"""
from ase.lattice import bulk
from ase.lattice.cubic import Diamond
from ase.constraints import FixAtoms
import ase.units as units

from quippy import set_fortran_indexing
from quippy.atoms import Atoms
from quippy.potential import Potential, Minim
from quippy.elasticity import youngs_modulus, poisson_ratio
from quippy.io import write
from quippy.crack import (print_crack_system,
                          G_to_strain,
                          thin_strip_displacement_y,
                          find_crack_tip_stress_field)

import atomeye

crack_direction = (-2, 1, 1)      # Miller index of x-axis
cleavage_plane = (1, 1, 1)        # Miller index of y-axis
crack_front = (0, 1, -1)          # Miller index of z-axis

width = 200.0*units.Ang              # Width of crack slab
height = 100.0*units.Ang             # Height of crack slab
vacuum = 100.0*units.Ang             # Amount of vacuum around slab
crack_seed_length = 40.0*units.Ang   # Length of seed crack
strain_ramp_length = 30.0*units.Ang  # Distance over which strain is ramped up
initial_G = 5.0*(units.J/units.m**2) # Initial energy flow to crack tip

relax_fmax = 0.025*units.eV/units.Ang  # Maximum force criteria for relaxation

output_file = 'crack.xyz'            # File to which structure will be written

set_fortran_indexing(False)

si_bulk = bulk('Si', 'diamond', a=5.431, cubic=True)
mm_pot = Potential("IP SW")
si_bulk.set_calculator(mm_pot)

unit_slab = Diamond(directions=[crack_direction,
                                cleavage_plane,
                                crack_front],
                    size=(1, 1, 1),
                    symbol='Si',
                    pbc=True,
                    latticeconstant=5.431)


unit_slab.positions[:, 1] += (unit_slab.positions[1, 1] -
                              unit_slab.positions[0, 1]) / 2.0
unit_slab.set_scaled_positions(unit_slab.get_scaled_positions())

surface = unit_slab.copy()
surface.center(vacuum, axis=1)

nx = int(width / unit_slab.cell[0, 0])
ny = int(height / unit_slab.cell[1, 1])

# make sure ny is even so slab is centered on a bond
if ny % 2 == 1:
    ny += 1

# make a supercell of unit_slab
crack_slab = unit_slab * (nx, ny, 1)

# open up the cell along x and y by introducing some vaccum
crack_slab.center(vacuum, axis=0)
crack_slab.center(vacuum, axis=1)



# centre the slab on the origin
crack_slab.positions[:, 0] -= crack_slab.positions[:, 0].mean()
crack_slab.positions[:, 1] -= crack_slab.positions[:, 1].mean()

orig_width = (crack_slab.positions[:, 0].max() -
              crack_slab.positions[:, 0].min())
orig_height = (crack_slab.positions[:, 1].max() -
               crack_slab.positions[:, 1].min())

atomeye.view(crack_slab)
while 1 ==1:
    pass

# print(('Made slab with %d atoms, original width and height: %.1f x %.1f A^2' %
#        (len(crack_slab), orig_width, orig_height)))
#
# top = crack_slab.positions[:, 1].max()
# bottom = crack_slab.positions[:, 1].min()
# left = crack_slab.positions[:, 0].min()
# right = crack_slab.positions[:, 0].max()
#
# # fix atoms in the top and bottom rows
# fixed_mask = ((abs(crack_slab.positions[:, 1] - top) < 1.0) |
#               (abs(crack_slab.positions[:, 1] - bottom) < 1.0))
# const = FixAtoms(mask=fixed_mask)
# crack_slab.set_constraint(const)
# print('Fixed %d atoms\n' % fixed_mask.sum())
#
#
# # ****** Apply initial strain ramp *****
#
# strain = G_to_strain(initial_G, E, nu, orig_height)
#
# crack_slab.positions[:, 1] += thin_strip_displacement_y(
#                                  crack_slab.positions[:, 0],
#                                  crack_slab.positions[:, 1],
#                                  strain,
#                                  left + crack_seed_length,
#                                  left + crack_seed_length + strain_ramp_length)
#
# print('Applied initial load: strain=%.4f, G=%.2f J/m^2' %
#       (strain, initial_G / (units.J / units.m**2)))
#
#
# # ***** Relaxation of crack slab  *****
#
# # optionally, relax the slab, keeping top and bottom rows fixed
# print('Relaxing slab...')
# crack_slab.set_calculator(mm_pot)
# minim = Minim(crack_slab, relax_positions=True, relax_cell=False)
# minim.run(fmax=relax_fmax)
#
# # Find initial position of crack tip
# crack_pos = find_crack_tip_stress_field(crack_slab, calc=mm_pot)
# print 'Found crack tip at position %s' % crack_pos
#
# # Save all calculated materials properties inside the Atoms object
# crack_slab.info['nneightol'] = 1.3 # nearest neighbour tolerance
# crack_slab.info['LatticeConstant'] = a0
# crack_slab.info['C11'] = c[0, 0]
# crack_slab.info['C12'] = c[0, 1]
# crack_slab.info['C44'] = c[3, 3]
# crack_slab.info['YoungsModulus'] = E
# crack_slab.info['PoissonRatio_yx'] = nu
# crack_slab.info['SurfaceEnergy'] = gamma
# crack_slab.info['OrigWidth'] = orig_width
# crack_slab.info['OrigHeight'] = orig_height
# crack_slab.info['CrackDirection'] = crack_direction
# crack_slab.info['CleavagePlane'] = cleavage_plane
# crack_slab.info['CrackFront'] = crack_front
# crack_slab.info['strain'] = strain
# crack_slab.info['G'] = initial_G
# crack_slab.info['CrackPos'] = crack_pos
# crack_slab.info['is_cracked'] = False
#
#
# # ******** Save output file **********
#
# # save results in extended XYZ format, including extra properties and info
# print('Writing crack slab to file %s' % output_file)
# write(output_file, crack_slab)
