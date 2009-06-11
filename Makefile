DYNARE_VERSION=4.0.4

DYNAREBASE=dynare-$(DYNARE_VERSION)

all:
.PHONY: all

check:
	make -C tests $@ DYNARE_VERSION=$(DYNARE_VERSION)
.PHONY: check

srctarball:
	make -C preprocessor clean
	make -C doc clean
	make -C tests clean
	rm -f matlab/dynare_m matlab/dynare_m.exe
	rm -f mex/2007a/* mex/2007b/* mex/octave/*.mex
	rm -f windows/*.exe
	tar cvzf ../$(DYNAREBASE).tar.gz --transform 's,^\./,$(DYNAREBASE)/,' --exclude=debian --exclude='*~' --exclude-vcs .
.PHONY: srctarball
