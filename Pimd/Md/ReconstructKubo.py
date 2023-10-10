import os
from pathlib import Path

import numpy as np
from ase.io import read, Trajectory
from ase.units import kB
from ase.calculators.lammpsrun import LAMMPS
from ase.calculators.singlepoint import SinglePointCalculator


root = Path().absolute()
temperature = 38
nbeads = 64
os.environ["ASE_LAMMPSRUN_COMMAND"] = "lmp_mpi"

azizfile = (root / "HeAziz.table").as_posix()
pair_style = "table spline 1000"
pair_coeff = [f"* * {azizfile} HeHe_Aziz"]
calc = LAMMPS(pair_style=pair_style, pair_coeff=pair_coeff)

print("Getting all traj in a list")
allconfs = []
for i in range(1, nbeads+1):
    trajfile = root / f"Traj/mlmd.traj_{i:02d}"
    confs = read(trajfile, index=":")
    allconfs.append(confs)
nsteps = len(allconfs[0])
print("Done")


def compute_centroid_atoms(confs, temperature):
    """
    Function to compute the centroid

    Parameters
    ----------

    confs: :class:`list` of :class:`ase.Atoms`
        The configurations of the quantum polymer
    temperature: :class:`float`

    Returns
    -------

    atoms: :class:`ase.Atoms`
        The centroid of the quantum polymer
    """
    nbead = len(confs)
    atoms = confs[0].copy()
    natoms = len(atoms)
    masses = confs[0].get_masses()
    momenta = np.zeros((natoms, 3))

    kBT = kB * temperature

    pos = np.zeros((nbead, natoms, 3))
    forces = np.zeros((nbead, natoms, 3))
    epot = np.zeros(nbead)
    stress = np.zeros((nbead, 6))
    cell = np.zeros((nbead, 3, 3))
    for ibead, at in enumerate(confs):
        pos[ibead] = at.get_positions()
        forces[ibead] = at.get_forces()
        epot[ibead] = at.get_potential_energy()
        stress[ibead] = at.get_stress()
        cell[ibead] = at.get_cell()

    cpos = pos.mean(axis=0)
    cforces = forces.mean(axis=0)
    cekin = 1.5 * natoms * kBT - 0.5 * np.sum((pos - cpos) * forces) / nbead
    cenergy = epot.mean()
    cstress = stress.mean(axis=0)
    ccell = cell.mean(axis=0)
    # We set all momenta to zero and put it on one atom in the z axis
    # This is just to be able to access the temperature through
    # atoms.get_temperature()
    momenta[0, 0] = np.sqrt(2*cekin*masses[0])

    atoms.set_positions(cpos)
    atoms.set_cell(ccell, True)
    atoms.set_momenta(momenta)

    calc = SinglePointCalculator(atoms,
                                 energy=cenergy,
                                 forces=cforces,
                                 stress=cstress)
    atoms.calc = calc
    return atoms


print("Computing Kubo")
kubo = Trajectory("TrajKubo.traj", "w")
for i in range(nsteps):
    if i % 200 == 0:
        print(f"Step {i+1}")
    beadconfs = []
    for j in range(nbeads):
        at = allconfs[j][i]
        at.calc = calc
        at.get_potential_energy()
        beadconfs.append(at)
    kuboat = compute_centroid_atoms(beadconfs, temperature)
    kubo.write(kuboat)
print("Done")

print(len(allconfs))
print(len(allconfs[0]))
