# Nelder-Mead optimization for BTE inversion (2D-Dots)
The single Al_Si_2D_dots_bkwrd_optim_NM.cpp file uses Nelder-Mead (NM) algorithm to find phonon properties using temperature measurements of the 2D-dots experiment
by solving the forward problem, the Boltzmann transport equation (BTE), using Monte Carlo (MC) simulations.

- material data files to read: BVK_Al.txt, BVK_Si.txt
- temperature measurement files to read: Al_Si_bkwrd**.txt
- the file requires boost package to generate random numbers. Add the package to the directory and run the cpp file with that package; use "g++ -I [boost_file] Al_Si_2D_dots_bkwrd_optim_NM.cpp -o [executable_file_name]"

Reference: 

M. Forghani et al., “Reconstruction of phonon relaxation times from systems featuring interfaces with unknown properties”, Physical Review B 97, 195440 (2018). https://journals.aps.org/prb/abstract/10.1103/PhysRevB.97.195440
