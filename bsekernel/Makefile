# Kernel for the BSE interpolation from BerkeleyGW
# Kernel originally by Felipe H. da Jornada (jornada@berkeley.edu)

# See common.mk for the actual make rules

.SUFFIXES:

DEFAULT all: mic host
clean: clean-mic clean-host

mic:
	mkdir -p mic
	cd mic && make -f ../common.mk
host:
	mkdir -p host
	cd host && make -f ../common.mk

clean-mic:
	cd mic && rm -f *.o *.mod *.optrpt *.x
clean-host:
	cd host && rm -f *.o *.mod *.optrpt *.x

.PHONY: mic host clean-mic clean-host

tarball:
	$(eval d := $(shell date +%Y%m%d))
	$(eval r := $(shell git rev-list --branches master master | wc -l))
	$(eval fname := "kernel_BSE-r${r}-${d}.tgz")
	echo "`git rev-list --max-count=1 --header master`" > VERSION
	tar --transform 's,^,kernel_BSE/,' -czf ${fname} *.f90 Makefile VERSION
	rm VERSION
