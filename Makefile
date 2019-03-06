ifeq ($(FC),gfortran)
	FFLAGS= -g -fdefault-real-8 -fdefault-double-8
	#FFLAGS= -g 
endif

ifeq ($(FC),ifort)
	FFLAGS= -g -fcheck=all
endif

OMP = -fopenmp

default: all

all: cspeed getoptf kd_tree runner

cspeed: cspeed.c
	gcc -c cspeed.c

getoptf: getoptf.f90
	$(FC) $(FFLAGS) -c getoptf.f90

kd_tree: kd_tree.f90
	$(FC) $(FFLAGS) -c kd_tree.f90

runner: runner.f90 kd_tree.o getoptf.o cspeed.o
	$(FC) $(FFLAGS) -o kd_tests runner.f90 getoptf.o kd_tree.o cspeed.o

clean:
	rm -rf *.mod *.o kd_tests
