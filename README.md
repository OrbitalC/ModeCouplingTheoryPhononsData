# Data for the article Mode-coupling theory of lattice dynamics for classical and quantum crystals


In this repository, I compiled the necessary files to recreate the results in the paper "Mode-coupling theory of lattice dynamics for classical and quantum crystals"

To run the simulations, you need 
- For the PIMD
    - At least 64 cpu available
    - A version of LAMMPS compiled with the replica PACKAGE
- For the spectral function
    - [https://github.com/tdep-developers/tdep](The TDEP code)
- For some pre/post-processing
    - Python with numpy, matplotlib, ase, h5py - can be installed with pip (not that for ASE, you need the developper version available on [https://gitlab.com/ase/ase](gitlab)

Note that to exactly reproduce the results, you need to slightly modify the sources of TDEP to compute the spectral function with equation 39 of the paper.
These modifications should happen in the file `src/lineshape/lineshape_helper.f90` and consist in modifying some omega with BigOmega(i).


For each approximation (Pimd, Classical, Scha, Perturbation), I've put the data in it's respective folder

For Pimd and Classical, you can find LAMMPS input to run the simulations in the Md/ folder.
To run them, you can simply use the bash script in increasing order. Note that you need to to have a binary lmp_mpi in your PATH, and in the case of PIMD, access to 64 CPU.
I also put the resulting trajectories in a ASE format.
In the Tdep folder, you will find everything to generate the spectral function on path in order to reproduce the figure 3. of the paper.
To do this, you can simply execute the bash files in increasing order. Nothe that you need to have added the tdep binary in your PATH, as proposed at the end of the TDEP compilation.

Once the spectral function have been computed for each approximation, you can execute the `PlotAllLineshape.py` script in the root of this repository to generate the figure 3 of the paper.


If you have any question, don't hesitate to raise an issue !
