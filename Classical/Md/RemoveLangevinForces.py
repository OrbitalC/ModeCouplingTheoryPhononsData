import os
from pathlib import Path

from ase.io import read, Trajectory
from ase.calculators.lammpsrun import LAMMPS


lmpbin = "lmp_mpi"
os.environ["ASE_LAMMPSRUN_COMMAND"] = lmpbin
root = Path().absolute()

aziz = root / "HeAziz.table"
pair_style = "table spline 1000"
pair_coeff = [f"* * {aziz} HeHe_Aziz"]

calc = LAMMPS(pair_style=pair_style, pair_coeff=pair_coeff)

# We need to do this to remove the Langevin forces to the results
oldtraj = read("mlmd.traj", index=":")
traj = Trajectory("MdTraj.traj", "w")
for at in oldtraj:
    at.calc = calc
    at.get_potential_energy()
    traj.write(at)
