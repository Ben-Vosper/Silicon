import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

x = [0.5, -0.5]
y = [0, 0]
z = [0, 0]

d_eq = abs(x[0]) + abs(x[1])
d_lat = d_eq*(2/3)*(2**0.5)
x_bb_rel = [d_lat/3, d_lat/3, d_lat/3,  -d_lat/3, -d_lat/3, -d_lat/3]
y_bb_rel = [d_lat, -d_lat, 0, -d_lat, d_lat, 0]
z_bb_rel = [-d_lat, -d_lat, d_lat, d_lat, d_lat, -d_lat]

x_bb = []
y_bb = []
z_bb = []
for q in range(6):
    if q < 3:
        x_bb.append(x_bb_rel[q] + x[0])
        y_bb.append(y_bb_rel[q] + y[0])
        z_bb.append(z_bb_rel[q] + z[0])
    else:
        x_bb.append(x_bb_rel[q] + x[1])
        y_bb.append(y_bb_rel[q] + y[1])
        z_bb.append(z_bb_rel[q] + z[1])

f = 1000

# Set up figure & 3D axis for animation
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1], projection='3d')
#ax.axis('off')
ax.set_xlim((-3, 3))
ax.set_ylim((-3, 3))
ax.set_zlim((-3, 3))

atoms = [ax.plot([], [], ls="none",  marker='o', markersize=20, color=(0, 0.5, 0.5), zorder=10),
         ax.plot([], [], ls="none",  marker='o', markersize=20, color=(0, 0.5, 0.5), zorder=10)]

bonds = [ax.plot([], [], ls="--", lw=2, color=(0, 1, 0.5))]

back_atoms = [ax.plot([], [], ls="none",  marker='o', markersize=20, color=(0.5,0.5, 0), zorder=9)]
back_bonds = [ax.plot([], [], ls="-", lw=2, color=(1, 0, 0.5))]
for q in range(5):
    back_atoms.append(ax.plot([], [], ls="none",  marker='o', markersize=20, color=(0.5,0.5, 0), zorder=9))
    back_bonds.append(ax.plot([], [], ls="-", lw=2, color=(1, 0, 0.5)))

# set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(30, 72)

# initialization function: plot the background of each frame
def init():
    elements = []
    for atom in atoms:
        atom = atom[0]
        atom.set_data([], [])
        atom.set_3d_properties([])
        elements.append(atom)

    for bond in bonds:
        bond = bond[0]
        bond.set_data([], [])
        bond.set_3d_properties([])
        elements.append(bond)

    for atom in back_atoms:
        atom = atom[0]
        atom.set_data([], [])
        atom.set_3d_properties([])
        elements.append(atom)

    for bond in back_bonds:
        bond = bond[0]
        bond.set_data([], [])
        bond.set_3d_properties([])
        elements.append(bond)

    return elements

# animation function.  This will be called sequentially with the frame number
def animate(i):
    elements = []
    x_offset = [-0.3*np.sin(2*np.pi*i/500), 0.3*np.sin(2*np.pi*i/500)]
    x_pos = np.add(x, x_offset)

    bondx = []
    bondy = []
    bondz = []
    for q in range(len(bonds)):
        bondx.append([x_pos[q], x_pos[q + 1]])
        bondy.append([y[q], y[q + 1]])
        bondz.append([z[q], z[q + 1]])
    for bond in bonds:
        q = bonds.index(bond)
        bond = bond[0]
        bond.set_data(bondx[q], bondy[q])
        bond.set_3d_properties(bondz[q])
        elements.append(bond)

    bondx = []
    bondy = []
    bondz = []
    for q in range(len(back_bonds)):
        if q < 3:
            bondx.append([x_bb[q], x_pos[0]])
            bondy.append([y_bb[q], y[0]])
            bondz.append([z_bb[q], z[0]])
        else:
            bondx.append([x_bb[q], x_pos[1]])
            bondy.append([y_bb[q], y[1]])
            bondz.append([z_bb[q], z[1]])
    for bond in back_bonds:
        q = back_bonds.index(bond)
        bond = bond[0]
        bond.set_data(bondx[q], bondy[q])
        bond.set_3d_properties(bondz[q])
        elements.append(bond)

    for atom in atoms:
        q = atoms.index(atom)
        atom = atom[0]
        atom.set_data([x_pos[q]], [y[q]])
        atom.set_3d_properties([z[q]])
        elements.append(atom)

    for atom in back_atoms:
        q = back_atoms.index(atom)
        atom = atom[0]
        atom.set_data([x_bb[q]], [y_bb[q]])
        atom.set_3d_properties([z_bb[q]])
        elements.append(atom)

    #ax.view_init(30, 72 + (0.1 * i))
    fig.canvas.draw()
    return elements

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=f, interval=1, blit=True)
plt.show()











