F90 = ifort -warn all -g -traceback -check all
LINK = ifort -g -traceback -check all
LIBS =

kernel_BSE.x: kernel_BSE.o common.o intkernel.o cells.o
	$(LINK) -o $@ $^ $(LIBS)

kernel_BSE.o: common.o intkernel.o
intkernel.o: common.o cells.o blas.o

.PHONY: clean
clean:
	rm -rf *.o *.mod *.x

%.o: %.f90
	$(F90) -c $<
