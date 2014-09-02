F90 = ftn
LINK = ftn
LIBS =

kernel_BSE.x: kernel_BSE.o common.o
	$(LINK) -o $@ $^ $(LIBS)

kernel_BSE.o: common.o intkernel.o
intkernel.o: common.o cells.o blas.o

.PHONY: clean
clean:
	rm -rf *.o *.mod *.x

%.o: %.f90
	$(F90) -c $<
