FC=`which gfortran`
FLAGS=#-check all -traceback -fpe0
SUBDIR=sub_prg/src

a.out : coulomb.o eigen.o angular_coup.o \
  	  lin_coup.o print_array.o special_fct.o \
	  global_constant.o input_data.o potentials.o \
        relative_potential.o coupling_matrix.o coupled_channels.o calc_profile.o \
	  scat_ccfull.o
	${FC} -o a.out coulomb.o eigen.o angular_coup.o \
	  lin_coup.o special_fct.o print_array.o \
	  global_constant.o input_data.o potentials.o \
        relative_potential.o coupling_matrix.o coupled_channels.o calc_profile.o \
	  scat_ccfull.o ${FLAGS}
coulomb.o : ${SUBDIR}/coulomb.f90
	${FC} ${FLAGS} -c ${SUBDIR}/coulomb.f90
eigen.o : ${SUBDIR}/eigen.f90
	${FC} ${FLAGS} -c ${SUBDIR}/eigen.f90
angular_coup.o : ${SUBDIR}/angular_coup.f90
	${FC} ${FLAGS} -c ${SUBDIR}/angular_coup.f90
lin_coup.o : ${SUBDIR}/lin_coup.f90
	${FC} ${FLAGS} -c ${SUBDIR}/lin_coup.f90
print_array.o : ${SUBDIR}/print_array.f90
	${FC} ${FLAGS} -c ${SUBDIR}/print_array.f90
special_fct.o : ${SUBDIR}/special_fct.f90
	${FC} ${FLAGS} -c ${SUBDIR}/special_fct.f90
global_constant.o : global_constant.f90
	${FC} ${FLAGS} -c global_constant.f90
input_data.o : input_data.f90
	${FC} ${FLAGS} -c input_data.f90
potentials.o : potentials.f90 input_data.o
	${FC} ${FLAGS} -c potentials.f90
relative_potential.o : relative_potential.f90 input_data.o potentials.o
	${FC} ${FLAGS} -c relative_potential.f90
coupling_matrix.o : coupling_matrix.f90 eigen.o angular_coup.o \
  				global_constant.o input_data.o \
				potentials.o eigen.o relative_potential.o
	${FC} ${FLAGS} -c coupling_matrix.f90
coupled_channels.o : coupled_channels.f90 lin_coup.o coulomb.o special_fct.o \
				global_constant.o input_data.o \
				relative_potential.o coupling_matrix.o
	${FC} ${FLAGS} -c coupled_channels.f90
calc_profile.o : calc_profile.f90 input_data.o potentials.o coupling_matrix.o \
  			relative_potential.o
	${FC} ${FLAGS} -c calc_profile.f90
scat_ccfull.o : scat_ccfull.f90 global_constant.o input_data.o \
			relative_potential.o coupling_matrix.o \
			coupled_channels.o calc_profile.o
	${FC} ${FLAGS} -c scat_ccfull.f90

.PHONY : clean
clean :
	$(RM) *.o *.mod


