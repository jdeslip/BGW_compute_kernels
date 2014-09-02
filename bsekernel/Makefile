# Kernel for the BSE interpolation from BerkeleyGW
# Kernel originally by Felipe H. da Jornada (jornada@berkeley.edu)

#F90 = ifort -warn all -g -traceback -check all -openmp
#LINK = ifort -g -traceback -check all -openmp

#FHJ: Note: I'm keeping -ip out of the compile flags.
#It takes forever to compile BGW with -ip!
F90 = ifort -warn all -openmp -O3 -no-prec-div -static -xHost

LINK = ifort -openmp
LIBS =

kernel_BSE.x: kernel_BSE.o common.o intkernel.o cells.o timing.o
	$(LINK) -o $@ $^ $(LIBS)

kernel_BSE.o: common.o intkernel.o timing.o
intkernel.o: common.o cells.o blas.o timing.o
timing.o: common.o

.PHONY: clean tarball
clean:
	rm -rf *.o *.mod *.x
tarball:
	tar --transform 's,^,kernel_BSE/,' -czf kernel_BSE.tgz *.f90 Makefile


%.o: %.f90
	$(F90) -c $<
