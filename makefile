FC=`which gfortran`
SUBDIR=./utils
OBJDIR=./obj
MODDIR=./modules

UTIL_LIB=${SUBDIR}/lib/libmylib.a
MYLIB=-I ${SUBDIR}/modules -L${SUBDIR}/lib -lmylib
FLAGS=-O2 -J $(MODDIR) ${MYLIB}

OBJ=${OBJDIR}/global_constant.o \
    ${OBJDIR}/input_data.o \
    ${OBJDIR}/potentials.o \
    ${OBJDIR}/relative_potential.o \
    ${OBJDIR}/coupling_matrix.o \
    ${OBJDIR}/coupled_channels.o \
    ${OBJDIR}/calc_profile.o \
    ${OBJDIR}/scat_ccfull.o

a.out: ${UTIL_LIB} ${OBJ}
	${FC} -o a.out ${OBJ} ${FLAGS}
${UTIL_LIB}: utils
	cd utils && make
${OBJDIR}/global_constant.o : global_constant.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/global_constant.o global_constant.f90
${OBJDIR}/input_data.o : input_data.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/input_data.o input_data.f90
${OBJDIR}/potentials.o : potentials.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/potentials.o potentials.f90
${OBJDIR}/relative_potential.o : relative_potential.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/relative_potential.o relative_potential.f90
${OBJDIR}/coupling_matrix.o : coupling_matrix.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/coupling_matrix.o coupling_matrix.f90
${OBJDIR}/coupled_channels.o : coupled_channels.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/coupled_channels.o coupled_channels.f90
${OBJDIR}/calc_profile.o : calc_profile.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/calc_profile.o calc_profile.f90
${OBJDIR}/scat_ccfull.o : scat_ccfull.f90
	${FC} ${FLAGS} -c -o ${OBJDIR}/scat_ccfull.o scat_ccfull.f90

.PHONY : clean
clean :
	$(RM) ${OBJDIR}/*.o ${MODDIR}/*.mod *~
