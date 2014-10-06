# Kernel for the Chi Summation from BerkeleyGW
# Borrowed from FHJ's BSE kernel.

#F90 = ifort -warn all -g -traceback -check all -openmp
#LINK = ifort -g -traceback -check all -openmp

F90 = ftn -openmp -fast -no-ipo #ifort -warn all -openmp -O3 -no-prec-div -static -xHost

LINK = ftn -fast -no-ipo -openmp # ifort -openmp
#LAPACKLIB = -L/opt/cray/libsci/12.1.01/PGI/121/mc12/lib/ -lsci_pgi_mp
LAPACKLIB="-mkl=parallel"
#LAPACKLIB=$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
#               $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a
LIBS = #$(LAPACKLIB)

kernel_chi_sum.x: chi_summation.o blas.o timing.o
	$(LINK) -o $@ $^ $(LAPACKLIB)

chi_summation.o: blas.o timing.o

.PHONY: clean tarball
clean:
	rm -rf *.o *.mod *.x
tarball:
	tar --transform 's,^,kernel_chi_sum/,' -czf kernel_chi_sum.tgz *.f90 Makefile


%.o: %.f90
	$(F90) -c $<
