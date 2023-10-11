from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ase.units import Hartree, Bohr


root = Path().absolute()
figname = "IrreducibleIFC.pdf"
# We will need this to convert tdep outpus
units = Hartree / Bohr**2

# Define some path to the results
pimdfd = root / "Pimd/Tdep"
clasfd = root / "Classical/Tdep"
pertfd = root / "Perturbation/Tdep"
schafd = root / "Scha/Tdep"

idx_max = 12  # 3
file = "outfile.irrifc_secondorder"
pimd = np.loadtxt(pimdfd / file)[:idx_max] * units
scha = np.loadtxt(schafd / file)[:idx_max] * units
clas = np.loadtxt(clasfd / file)[:idx_max] * units
harm = np.loadtxt(pertfd / file)[:idx_max] * units

idx_max = 5  # 5
file = "outfile.irrifc_thirdorder"
pimd3 = np.loadtxt(pimdfd / file)[:idx_max] * units / Bohr
scha3 = np.loadtxt(schafd / file)[:idx_max] * units / Bohr
clas3 = np.loadtxt(clasfd / file)[:idx_max] * units / Bohr
harm3 = np.loadtxt(pertfd / file)[:idx_max] * units / Bohr

# A little bit of ordering to have a nicer plot
idx = np.abs(pimd) > 1e-3
pimd = pimd[idx]
scha = scha[idx]
clas = clas[idx]
harm = harm[idx]
idx = np.argsort(pimd)
pimd = pimd[idx]
scha = scha[idx]
clas = clas[idx]
harm = harm[idx]

# The same here
idx3 = np.abs(pimd3) > 1e-3
pimd3 = pimd3[idx3]
scha3 = scha3[idx3]
clas3 = clas3[idx3]
harm3 = harm3[idx3]
idx3 = np.argsort(pimd3)
pimd3 = pimd3[idx3]
scha3 = scha3[idx3]
clas3 = clas3[idx3]
harm3 = harm3[idx3]


# And now the plots
fig = plt.figure(figsize=(20, 15), constrained_layout=True)
gd = GridSpec(2, 1, fig)

# Just to have something looking nice
mpl.rcParams["lines.linewidth"] = 5
mpl.rcParams["lines.markeredgecolor"] = "k"
mpl.rcParams["lines.markersize"] = 25
mpl.rcParams["lines.markeredgewidth"] = 5

mpl.rcParams["font.size"] = 30

mpl.rcParams["axes.linewidth"] = 5

mpl.rcParams["xtick.top"] = True
mpl.rcParams["xtick.major.size"] = 12
mpl.rcParams["xtick.major.width"] = 5
mpl.rcParams["xtick.direction"] = "in"

mpl.rcParams["ytick.right"] = True
mpl.rcParams["ytick.major.size"] = 12
mpl.rcParams["ytick.major.width"] = 5
mpl.rcParams["ytick.direction"] = "in"

ax0 = fig.add_subplot(gd[0])
ax1 = fig.add_subplot(gd[1])

ax0.plot(pimd, marker="o", label="Mode-coupling PIMD")
ax0.plot(scha, marker="o", label="SCHA")
ax0.plot(clas, marker="o", label="Mode-coupling MD")
ax0.plot(harm, marker="o", label="Perturbation theory")
ax0.text(0.025, 0.90, "a) Frequency matrix", transform=ax0.transAxes,
         color="k")

ax1.plot(pimd3, marker="o", label="Mode-coupling PIMD")
ax1.plot(scha3, marker="o", label="SCHA")
ax1.plot(clas3, marker="o", label="Mode-coupling MD")
ax1.plot(harm3, marker="o", label="Harmonic")
ax1.text(0.025, 0.90, r"b) 2$\mathrm{nd}$ order vertex",
         transform=ax1.transAxes, color="k")

ax0.set_xticks(np.arange(len(pimd)))
ax0.set_xlabel("Index")
ax0.set_ylabel(r"Coefficient [eV/$\mathring{A}^2$]")

ax1.set_xticks(np.arange(len(pimd3)))
ax1.set_xlabel("Index")
ax1.set_ylabel(r"Coefficient [eV/$\mathring{A}^3$]")

ax0.legend(ncol=2, loc=(0.025, 0.64))
fig.savefig(figname)
plt.show()
