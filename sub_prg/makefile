# This makefile produces a library file from programs in the directory "src/".
SRCDIR=src
MODDIR=module
OBJDIR=obj

CMPL=/usr/bin/ifort
FLAG=-fPIC -module $(MODDIR)

SHARED_LIB=libmylib.so

OBJ=$(OBJDIR)/coulomb.o \
	$(OBJDIR)/eigen.o \
	$(OBJDIR)/angular_coup.o \
	$(OBJDIR)/gaus_int.o \
	$(OBJDIR)/gaus_weights.o \
	$(OBJDIR)/lin_coup.o \
	$(OBJDIR)/print_array.o \
	$(OBJDIR)/special_fct.o \


lib : $(SHARED_LIB)

$(SHARED_LIB) : $(OBJ) 
	$(CMPL) -shared -o $(SHARED_LIB) $(OBJ)


$(OBJDIR)/coulomb.o : $(SRCDIR)/coulomb.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/coulomb.o \
	$(SRCDIR)/coulomb.f90
$(OBJDIR)/eigen.o : $(SRCDIR)/eigen.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/eigen.o \
	$(SRCDIR)/eigen.f90
$(OBJDIR)/angular_coup.o : $(SRCDIR)/angular_coup.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/angular_coup.o \
	$(SRCDIR)/angular_coup.f90
$(OBJDIR)/gaus_int.o : $(SRCDIR)/gaus_int.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/gaus_int.o \
	$(SRCDIR)/gaus_int.f90
$(OBJDIR)/gaus_weights.o : $(SRCDIR)/gaus_weights.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/gaus_weights.o \
	$(SRCDIR)/gaus_weights.f90
$(OBJDIR)/lin_coup.o : $(SRCDIR)/lin_coup.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/lin_coup.o \
	$(SRCDIR)/lin_coup.f90
$(OBJDIR)/print_array.o : $(SRCDIR)/print_array.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/print_array.o \
	$(SRCDIR)/print_array.f90
$(OBJDIR)/special_fct.o : $(SRCDIR)/special_fct.f90
	$(CMPL) -c $(FLAG) \
	-o $(OBJDIR)/special_fct.o \
	$(SRCDIR)/special_fct.f90

.PHONY : clean
clean :
	rm $(OBJ) $(MODDIR)/*
