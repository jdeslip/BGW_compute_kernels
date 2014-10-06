# Kernel for the Chi Summation from BerkeleyGW
# Borrowed from FHJ's BSE kernel.

#F90 = ifort -warn all -g -traceback -check all -openmp
#LINK = ifort -g -traceback -check all -openmp

F90 = ftn #ifort -warn all -openmp -O3 -no-prec-div -static -xHost

LINK = ftn # ifort -openmp
LAPACKLIB = -L/opt/cray/libsci/12.1.01/PGI/121/mc12/lib/ -lsci_pgi_mp
LIBS = #$(LAPACKLIB)

kernel_chi_sum.x: chi_summation.o blas.o timing.o
	$(LINK) -o $@ $^ $(LIBS)

chi_summation.o: blas.o timing.o

.PHONY: clean tarball
clean:
	rm -rf *.o *.mod *.x
tarball:
	tar --transform 's,^,kernel_chi_sum/,' -czf kernel_chi_sum.tgz *.f90 Makefile


%.o: %.f90
	$(F90) -c $<
