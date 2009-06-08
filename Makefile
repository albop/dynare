# Makefile for creating the source tarball

DYNAREBASE=dynare-4.0.4

srctarball:
	make -C preprocessor clean
	make -C doc clean
	rm -f matlab/dynare_m matlab/dynare_m.exe
	rm -f mex/2007a/* mex/2007b/* mex/octave/*.mex
	rm -f windows/*.exe
	tar cvzf ../$(DYNAREBASE).tar.gz --transform 's,^\./,$(DYNAREBASE)/,' --exclude=debian --exclude='*~' --exclude-vcs .
