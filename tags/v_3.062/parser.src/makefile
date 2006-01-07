all: dynare_m.exe dynare_s.exe dynare_g.exe dynare_OxPro.exe dynare_Oxcons.exe

precompiler.o: precompiler.cc precompiler.hh d_ll.c d_tab.c
	g++ -mno-cygwin -c precompiler.cc

dynare_mb: dynare_m.o
	gcc -mno-cygwin -g -o dynare_m.exe dynare_m.o d_tab.o d_ll.o -lefence
dynare_m.exe: dynare_m.o d_tab.o d_ll.o precompiler.o
	gcc -mno-cygwin -g -o dynare_m.exe dynare_m.o d_tab.o d_ll.o precompiler.o -lstdc++
	cp dynare_m.exe ../matlab

dynare_m.o: dynare.c d.h d_tab.h
	gcc -mno-cygwin -g  -c dynare.c -o dynare_m.o -DMATLAB -Wall #-DDEBUG

dynare_s.exe: dynare_s.o d_tab.o d_ll.o
	gcc -mno-cygwin -g  -o dynare_s.exe dynare_s.o d_tab.o d_ll.o
	cp dynare_s.exe ../scilab

dynare_s.o: dynare.c d.h d_tab.h
	gcc -mno-cygwin -g  -c dynare.c -o dynare_s.o -DSCILAB -Wall

dynare_g.exe: dynare_g.o d_tab.o d_ll.o
	gcc -mno-cygwin -g  -o dynare_g.exe dynare_g.o d_tab.o d_ll.o
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
	gcc -mno-cygwin -g -c d_ll.c


d_tab.o: d_tab.c
	gcc -mno-cygwin -g -c -DYYDEBUG=1 d_tab.c

d_tab.c d_tab.h: d.y d.h
	bison -d --verbose -od_tab.c d.y

d_ll.c: dyn.l d_tab.h
	flex -i -od_ll.c dyn.l


