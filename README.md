# Data for the article Mode-coupling theory of lattice dynamics for classical and quantum crystals


In this repository, I compiled the necessary files to recreate the results in the paper "Mode-coupling theory of lattice dynamics for classical and quantum crystals"

To run the simulations, you need 
- At least 64 cpu available
- A version of LAMMPS compiled with the replica PACKAGE
- The TDEP code
- Python with numpy, matplotlib, ase, h5py

Note that to exactly reproduce the results, you need to slightly modify the sources of TDEP to compute the spectral function with equation 39 of the paper.


For each approximation (Pimd, Classical, Scha, Perturbation), I've put the data in it's respective folder

For Pimd and Classical, you can find LAMMPS input to run the simulations in the Md/ folder.
You will also find a Tdep Folder with a python script to generate input files for Tdep and some bash script to run TDEP and extract the spectral function on a path.
Finally, the TdepReference folder presents the force constants and spectral function that were use to generate figure 3. of the paper.

I've also added the trajectory I used in the paper as ASE trajectory files.


If you use the script `PlotAllLineshape.py` in the root of this repository, you should reproduce the figure 3 of the paper.


If you have any question, don't hesitate to raise an issue !
