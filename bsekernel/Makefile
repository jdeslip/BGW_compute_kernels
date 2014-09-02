F90 = ftn
LINK = ftn
LIBS =

kernel_BSE.x: kernel_BSE.o common.o
	$(LINK) -o $@ $^ $(LIBS)

kernel_BSE.o: common.o intkernel.o
intkernel.o: common.o cells.o

%.o: %.f90
	$(F90) -c -o $@ $^
