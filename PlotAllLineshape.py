from pathlib import Path

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from h5py import File
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec


root = Path().absolute()
savefile = "Fcc_38K_Aziz.pdf"
ommax = 22
thz = 4.135665538536

# Load data from pimd
folder = root / "Pimd/Tdep"
ls = File(folder / "outfile.phonon_spectral_function.hdf5")
omega_pimd = np.array(ls.get("energy_values")) * thz
specf_pimd = np.array(ls.get("spectral_function")) / thz
qval_pimd = np.array(ls.get("q_values"))
lticks = [r"$\Gamma$", "X", "U/K", r"$\Gamma$", "L"]
qticks_pimd = np.array(ls.get("q_ticks"))

# Load data from classical md
folder = root / "Classical/Tdep"
omega_md = np.array(ls.get("energy_values")) * thz
specf_md = np.array(ls.get("spectral_function")) / thz
qval_md = np.array(ls.get("q_values"))
lticks = [r"$\Gamma$", "X", "U/K", r"$\Gamma$", "L"]
qticks_md = np.array(ls.get("q_ticks"))

# Load data from scha
folder = root / "Scha/Tdep"
ls = File(folder / "outfile.phonon_spectral_function.hdf5")
omega_sc = np.array(ls.get("energy_values")) * thz
specf_sc = np.array(ls.get("spectral_function")) / thz
qval_sc = np.array(ls.get("q_values"))
lticks = [r"$\Gamma$", "X", "U/K", r"$\Gamma$", "L"]
qticks_sc = np.array(ls.get("q_ticks"))

# Load data from perturbation theory
folder = root / "Perturbation/Tdep"
ls = File(folder / "outfile.phonon_spectral_function.hdf5")
omega_pt = np.array(ls.get("energy_values")) * thz
specf_pt = np.array(ls.get("spectral_function")) / thz
qval_pt = np.array(ls.get("q_values"))
lticks = [r"$\Gamma$", "X", "U/K", r"$\Gamma$", "L"]
qticks_pt = np.array(ls.get("q_ticks"))

# We add a little background to make a nicer plot
# It doesn't change the result, but allows to use logspace for the color
# Instead of something more convoluted
specf_pimd += 1e-2 / thz
specf_md += 1e-2 / thz
specf_sc += 1e-2 / thz
specf_pt += 1e-2 / thz


# ========================================================== #
# Load experimental data

t_100 = np.loadtxt(root / "Expt/T_100.dat")
l_100 = np.loadtxt(root / "Expt/L_100.dat")
t1_110 = np.loadtxt(root / "Expt/T1_110.dat")
t2_110 = np.loadtxt(root / "Expt/T2_110.dat")
l_110 = np.loadtxt(root / "Expt/L_110.dat")
t_111 = np.loadtxt(root / "Expt/T_111.dat")
l_111 = np.loadtxt(root / "Expt/L_111.dat")

l_100[:, 0] *= qval_pimd[100]
t_100[:, 0] *= qval_pimd[100]

t_111[:, 0] *= 2.0 * (qval_pimd[-1] - qval_pimd[300])
t_111[:, 0] += qval_pimd[300]
l_111[:, 0] *= 2.0 * (qval_pimd[-1] - qval_pimd[300])
l_111[:, 0] += qval_pimd[300]

t1_110[:, 0] *= qval_pimd[100] - qval_pimd[300]
t1_110[:, 0] += qval_pimd[300]
t2_110[:, 0] *= qval_pimd[100] - qval_pimd[300]
t2_110[:, 0] += qval_pimd[300]
l_110[:, 0] *= qval_pimd[100] - qval_pimd[300]
l_110[:, 0] += qval_pimd[300]

# ========================================================== #

fig = plt.figure(figsize=(40, 30), constrained_layout=True)
gd = GridSpec(2, 31, fig)

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


ax0 = fig.add_subplot(gd[0, :15])
ax1 = fig.add_subplot(gd[0, 15:30], sharex=ax0, sharey=ax0)
ax2 = fig.add_subplot(gd[1, :15], sharex=ax0, sharey=ax0)
ax3 = fig.add_subplot(gd[1, 15:30], sharex=ax0, sharey=ax0)
ax4 = fig.add_subplot(gd[:, -1])

gx, gy = np.meshgrid(qval_pimd, omega_pimd)

cmap = "RdYlBu_r"  # sns.color_palette("light:b", as_cmap=True)
green = "#2ca02c"

norm = LogNorm(vmin=specf_pimd.min(), vmax=specf_pimd.max())
ax0.pcolormesh(gx, gy, specf_pimd, norm=norm, cmap=cmap,
               shading="gouraud", rasterized=True)

ax0.plot(t_100[:, 0], t_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax0.plot(l_100[:, 0], l_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax0.plot(t_111[:, 0], t_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax0.plot(l_111[:, 0], l_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax0.plot(t1_110[:, 0], t1_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax0.plot(t2_110[:, 0], t2_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax0.plot(l_110[:, 0], l_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax0.text(0.025, 0.95, "a) Mode-coupling quantum", transform=ax0.transAxes,
         color="w")


gx, gy = np.meshgrid(qval_md, omega_md)

cmap = "RdYlBu_r"  # sns.color_palette("light:b", as_cmap=True)

cb = ax1.pcolormesh(gx, gy, specf_md, norm=norm, cmap=cmap,
                    shading="gouraud", rasterized=True)

ax1.plot(t_100[:, 0], t_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax1.plot(l_100[:, 0], l_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax1.plot(t_111[:, 0], t_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax1.plot(l_111[:, 0], l_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax1.plot(t1_110[:, 0], t1_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax1.plot(t2_110[:, 0], t2_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax1.plot(l_110[:, 0], l_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax1.text(0.025, 0.95, "b) Mode-coupling classical", transform=ax1.transAxes,
         color="w")


gx, gy = np.meshgrid(qval_sc, omega_sc)

cmap = "RdYlBu_r"  # sns.color_palette("light:b", as_cmap=True)

ax2.pcolormesh(gx, gy, specf_sc, norm=norm, cmap=cmap,
               shading="gouraud", rasterized=True)

ax2.plot(t_100[:, 0], t_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax2.plot(l_100[:, 0], l_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax2.plot(t_111[:, 0], t_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax2.plot(l_111[:, 0], l_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax2.plot(t1_110[:, 0], t1_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax2.plot(t2_110[:, 0], t2_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax2.plot(l_110[:, 0], l_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax2.text(0.025, 0.95, "c) SCHA", transform=ax2.transAxes,
         color="w")

gx, gy = np.meshgrid(qval_pt, omega_pt)

ax3.pcolormesh(gx, gy, specf_sc, norm=norm, cmap=cmap,
               shading="gouraud", rasterized=True)

ax3.plot(t_100[:, 0], t_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax3.plot(l_100[:, 0], l_100[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax3.plot(t_111[:, 0], t_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax3.plot(l_111[:, 0], l_111[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax3.plot(t1_110[:, 0], t1_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax3.plot(t2_110[:, 0], t2_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)
ax3.plot(l_110[:, 0], l_110[:, 1], ls="", marker="o",
         markersize=12.5, markeredgewidth=1.5, c=green)

ax3.text(0.025, 0.95, "d) Perturbation theory", transform=ax3.transAxes,
         color="w")


plt.colorbar(cb, ax4)
ax4.set_ylabel(r"$\chi''(\omega) [meV^{-1}]$")

ax0.set_xlim(qval_pimd.min(), qval_pimd.max())
ax0.set_ylim(omega_pimd.min(), omega_pimd.max())
ax0.set_ylim(omega_pimd.min(), ommax)
ax0.set_xticks(qticks_pimd)
ax0.set_xticklabels(lticks)
ax0.set_ylabel("Frequency [meV]")
ax1.set_ylabel("Frequency [meV]")
ax2.set_ylabel("Frequency [meV]")
ax3.set_ylabel("Frequency [meV]")

fig.savefig(savefile)
plt.show()
