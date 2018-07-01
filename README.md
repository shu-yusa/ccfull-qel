# ccfull_qel
Code for calculating quasi-elastic scattering using the coupled-channels method.

## Prerequisites
`make` and `gfortran` need to be installed.
For mac, you can install `gfortran` by Homebrew
```bash
brew install gcc
brew install gfortran
```
or you can install from official binary [GFortranBinaries - GCC Wiki](https://gcc.gnu.org/wiki/GFortranBinaries).


## Compile & Run
Compile the source codes by the following.
```bash
make
```
This will produce an executable file `a.out`. Run calculation by
```bash
./a.out
```

## Input file
`input_scat` is the input file for the calculation. You can configure reaction details and output directory.

## Details of source codes and modules
* `global_constant.f90` ... Definition of constants.
* `input_data.f90` ... Read input parameters from the input file.
* `scat_ccfull.f90` ... The main program.
* `potentials.f90` ... Potential class is defined. Methods for calculating Coulomb, nuclear, and centrifugal potentials are defined.
* `relative_potential.f90` ... Class for calculating the relative potential for the reaction.
* `coupling_matrix.f90` ... Class for calculating the coupling matrix for the reaction.
* `coupled_channels.f90` ... Class for calculating coupled-channels calculation by Numerov method. 
* `calc_profile.f90` ... Class for output calculation details into a file `calc_info`.
* `utils/` ... Utility functions are contained in this directory.

## Output
Calculation results are output in the directory specified in the input file. 
Outputs are consisted of the following files.
* `calc_info` ... Calculation configurations.
* `fusion.dat` ... Fusion cross sections.
* `fus_bar_dist.dat` ... Fusion barrier distributions.
* `pot_nucl.dat` ... Nuclear potentions.
* `potential.dat` ... Real and imaginary parts of the total potential.
* `angular/` ... Quasi-elastic scattering cross sections and barrier distributions at specified angles.
* `angular_dist/` ... Angular distribution of elastic differential scattering cross sections.
* `Q_val_dist/` ... Q-value distribution of the scattering.

## References
[K. Hagino, N. Rowley, and A.T. Kruppa, Comput. Phys. Comm 123 (1999) 143 - 152](http://www.nucl.phys.tohoku.ac.jp/~hagino/ccfull.pdf)  
[Tohoku University Repository（TOUR）](https://tohoku.repo.nii.ac.jp/?action=pages_view_main&active_action=repository_view_main_item_detail&item_id=70405&item_no=1&page_id=33&block_id=38)
