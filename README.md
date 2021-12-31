# Nelder-Mead optimization for BTE inversion (1D-TTG)
The single Al_Si_2D_dots_bkwrd_optim_NM.cpp file uses Nelder-Mead (NM) algorithm to find phonon properties using temperature measurements of the 2D-dots experiment
by solving the forward problem, the Boltzmann transport equation (BTE), using Monte Carlo (MC) simulations.

- material data files to read: BVK_Al.txt, BVK_Si.txt
- temperature measurement files to read: Si_cos_bkwrd**.txt
- the file requires boost package to generate random numbers. Add the package to the directory and run the cpp file with that package; use "g++ -I [boost_file] Si_cos_bkwrd_optim_NM.cpp -o [executable_file_name]"