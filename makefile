MUMPS = $(HOME)/pak/libs-mcm-large/MUMPS_5.1.2
include $(MUMPS)/Makefile.inc

CC=gcc-4.4
CPP=g++-4.4
FC=gfortran-4.4

MPICH = $(HOME)/pak/libs-mcm-large/mpich-3.2.1/lib/.libs
MPICHINC = $(MPICH)/../../src/include
MUMPSINC = $(MUMPS)/include

FLAGS= -fPIC -O0 #-fopenmp
INC=-I./src -I$(MUMPSINC) -I$(MPICHINC)
FLIB = $(FLIBPAK) $(MPICH)/libmpifort.a $(MPICH)/libmpi.a -L/usr/lib64/ -lpthread
FLIBPAK = $(MUMPS)/lib/libdmumps.a $(MUMPS)/lib/libmumps_common.a $(SCALAP) $(HOME)/pak/libs-mcm-large/lapack-3.7.1/librefblas.a $(HOME)/pak/libs-mcm-large/lapack-3.7.1/liblapack.a -L/usr/lib64/ -lpthread $(MUMPS)/lib/libpord.a

SRC=PAKS
ODIRF=obj
ODIRF90=obj
ODIRC=obj
ODIRCPP=obj

_OBJF90= $(patsubst %.f90,%.o,$(subst $(SRC)/,,$(wildcard $(SRC)/*.f90)))
OBJF90 = $(patsubst %,$(ODIRF90)/%,$(_OBJF90))

_OBJF:= $(patsubst %.for,%.o,$(subst $(SRC)/,,$(wildcard $(SRC)/*.for)))
OBJF = $(patsubst %,$(ODIRF)/%,$(_OBJF))

_OBJC= $(patsubst %.c,%.o,$(subst $(SRC)/,,$(wildcard $(SRC)/*.c)))
OBJC = $(patsubst %,$(ODIRC)/%,$(_OBJC))

_OBJCPP= $(patsubst %.cpp,%.o,$(subst $(SRC)/,,$(wildcard $(SRC)/*.cpp)))
OBJCPP = $(patsubst %,$(ODIRCPP)/%,$(_OBJCPP))

all: PakS

PakS: $(OBJF90) $(OBJF) $(OBJC) $(OBJCPP) makefile
	$(FC)  $(FLAGS) -static $(OBJF90) $(OBJF) $(OBJC) $(OBJCPP) $(FLIB) -lstdc++ -o paks
	
$(ODIRF90)/%.o: $(SRC)/%.f90
	$(FC) $(FLAGS) $(INC) -c -o $@ $<
	
$(ODIRF)/%.o: $(SRC)/%.for
	$(FC) $(FLAGS) $(INC) -c -o $@ $<
	
$(ODIRC)/%.o: $(SRC)/%.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	
$(ODIRCPP)/%.o: $(SRC)/%.cpp
	$(CPP) $(CFLAGS) $(INC) -c -o $@ $<
	
clean:
	rm obj/*.o
	rm *.mod
	rm paks

# DO NOT DELETE THIS LINE -- make depend depends on it.
