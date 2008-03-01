CC = gcc
CPP = g++

CPPFLAGS = -Wall

ifeq ($(shell uname -o), Cygwin)
	# Detection of uninitialized variables is buggy in Cygwin and generates spurious warnings
	CPPFLAGS += -Wno-uninitialized
	CPPFLAGS += -mno-cygwin
endif

ifeq ($(CROSS_WIN32), yes)
	CC = i586-mingw32msvc-cc
	CPP = i586-mingw32msvc-g++
	# Detection of uninitialized variables is buggy in MinGW and generates spurious warnings
	CPPFLAGS += -Wno-uninitialized
endif

ifeq ($(DEBUG),yes)
	CPPFLAGS += -ggdb
else
	CPPFLAGS += -O3
endif

ifeq ($(VALGRIND), yes)
	CPPFLAGS = -Wall -O -g -fno-inline
endif

all: dynare_m.exe dynare_s.exe

precompiler.o: precompiler.cc precompiler.hh d_ll.c d_tab.c
	$(CPP) $(CPPFLAGS) -c precompiler.cc

dynare_m.exe: dynare_m.o d_tab.o d_ll.o precompiler.o
	$(CC) $(CPPFLAGS) -o dynare_m.exe dynare_m.o d_tab.o d_ll.o precompiler.o -lstdc++
	cp dynare_m.exe ../matlab

dynare_m.o: dynare.c d.h d_tab.h
	$(CC) $(CPPFLAGS)  -c dynare.c -o dynare_m.o -DMATLAB

dynare_s.exe: dynare_s.o d_tab.o d_ll.o  precompiler.o
	$(CC) $(CPPFLAGS)  -o dynare_s.exe dynare_s.o d_tab.o d_ll.o precompiler.o -lstdc++
	cp dynare_s.exe ../scilab

dynare_s.o: dynare.c d.h d_tab.h
	$(CC) $(CPPFLAGS)  -c dynare.c -o dynare_s.o -DSCILAB

dynare_g.exe: dynare_g.o d_tab.o d_ll.o  precompiler.o
	gcc -mno-cygwin -g  -o dynare_g.exe dynare_g.o d_tab.o d_ll.o precompiler.o
	cp dynare_g.exe ../gauss

dynare_g.o: dynare.c d.h d_tab.h
	gcc -mno-cygwin -g  -c dynare.c -o dynare_g.o -DGAUSS -Wall #-DDEBUG

dynare_OxPro.exe: dynare_OxPro.o d_tab.o d_ll.o
	gcc -mno-cygwin -g  -o dynare_OxPro.exe dynare_OxPro.o d_tab.o d_ll.o
	cp dynare_OxPro.exe ../gauss

dynare_OxPro.o: dynare.c d.h d_tab.h
	gcc -mno-cygwin -g  -c dynare.c -o dynare_OxPro.o -DGAUSS -DOXGAUSS -Wall

dynare_Oxcons.exe: dynare_Oxcons.o d_tab.o d_ll.o
	gcc -mno-cygwin -g  -o dynare_Oxcons.exe dynare_Oxcons.o d_tab.o d_ll.o
	cp dynare_Oxcons.exe ../gauss

dynare_Oxcons.o: dynare.c d.h d_tab.h
	gcc -mno-cygwin -g  -c dynare.c -o dynare_Oxcons.o -DGAUSS -DOXGAUSS -DOXGAUSS_GNU -Wall

d_ll.o: d_ll.c d_tab.h
	$(CC) $(CPPFLAGS) -c d_ll.c


d_tab.o: d_tab.c
	$(CC) $(CPPFLAGS) -c -DYYDEBUG=1 d_tab.c

d_tab.c d_tab.h: d.y d.h
	bison -d --verbose -od_tab.c d.y

d_ll.c: dyn.l d_tab.h
	flex -i -od_ll.c dyn.l


