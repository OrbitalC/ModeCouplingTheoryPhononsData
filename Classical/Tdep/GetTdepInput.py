from pathlib import Path

import numpy as np
from ase.io import read
from ase.units import GPa


root = Path().absolute()
nthrow = 100
traj = read(root.parent / "Md/Traj.traj", f"{nthrow}:")
sim_temperature = 38
dt = 0.5  # in fs

nat = len(traj[0])
nsteps = len(traj)

pos = np.zeros((nsteps, nat, 3))
forces = np.zeros((nsteps, nat, 3))
statdata = np.zeros((nsteps, 13))
for i, at in enumerate(traj):
    pos[i] = at.get_scaled_positions()
    forces[i] = at.get_forces()
    stress = at.get_stress() / GPa

    statdata[i, 0] = i + 1
    statdata[i, 1] = i * dt
    statdata[i, 2] = at.get_total_energy()
    statdata[i, 3] = at.get_potential_energy()
    statdata[i, 4] = at.get_kinetic_energy()
    statdata[i, 5] = at.get_temperature()
    statdata[i, 6] = -stress[:3].sum() / 3.0
    statdata[i, 7:] = stress

fmt = "%20.15f"
np.savetxt("infile.positions", pos.reshape((nat * nsteps, 3)), fmt)
np.savetxt("infile.forces", forces.reshape((nat * nsteps, 3)), fmt)
fmt = "% 4d" + 12 * "%15.10f"
np.savetxt("infile.stat", statdata, fmt)

with open("infile.meta", "w") as fd:
    fd.write(f"{nat:6d}  # Number of atoms\n")
    fd.write(f"{nsteps:6d}  # Number of timestep\n")
    fd.write(f"{dt:6.2f}  # Timestep in fs, useless here\n")
    fd.write(f"{sim_temperature:6.2f}  # Temperature, in K")
