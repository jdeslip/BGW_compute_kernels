F90_MIC = mpiifort -warn all -openmp -O3 -no-prec-div -mmic -g -qopt-report=5 #-qopt-matmul #-qnolib-inline
F90_HOST = mpiifort -warn all -openmp -O3 -no-prec-div -static -xAVX -g -qopt-report=5

LINK_MIC = mpiifort -openmp -mmic
LINK_HOST = mpiifort -openmp -xAVX

LIBS_MIC =
LIBS_HOST =

ARCH = $(notdir $(CURDIR))
ifeq ($(ARCH),mic)
  F90 = $(F90_MIC)
  LINK = $(LINK_MIC)
else ifeq ($(ARCH),host)
  F90 = $(F90_HOST)
  LINK = $(LINK_HOST)
endif

kernel_BSE.$(ARCH).x: kernel_BSE.o common_m.mod intkernel_m.mod cells_m.mod timing_m.mod blas_m.mod
	$(LINK) -o $@ $(patsubst %_m.mod,%.o,$^) $(LIBS)

kernel_BSE.o: common_m.mod intkernel_m.mod timing_m.mod
intkernel_m.mod: common_m.mod cells_m.mod blas_m.mod timing_m.mod
cells_m.mod: common_m.mod
timing_m.mod: common_m.mod

%.o %_m.mod: ../%.f90
	$(F90) -c $<
