#F90 = ifort -warn all -g -traceback -check all -openmp
#LINK = ifort -g -traceback -check all -openmp
F90 = ifort -warn all -fast -no-ip -openmp
LINK = ifort -fast -no-ip -openmp
LIBS =

kernel_BSE.x: kernel_BSE.o common.o intkernel.o cells.o timing.o
	$(LINK) -o $@ $^ $(LIBS)

kernel_BSE.o: common.o intkernel.o timing.o
intkernel.o: common.o cells.o blas.o timing.o
timing.o: common.o

.PHONY: clean
clean:
	rm -rf *.o *.mod *.x

%.o: %.f90
	$(F90) -c $<
