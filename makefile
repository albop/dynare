DATE = $(shell date +%m%d)
VERSION = v3_065
DYNARE_ROOT = dynare_v3
FNAME_MAT = $(DYNARE_ROOT)/zips/dyn_mat_$(VERSION)
FNAME_PAR = $(DYNARE_ROOT)/zips/dyn_parser_$(VERSION)
MZIP_FILE = dyn_mat_$(VERSION).zip 
PZIP_FILE = dyn_parser_$(VERSION).zip 

all: install

parser.src/dynare_m.exe: parser.src/dynare.c parser.src/d.y parser.src/dyn.l doc/guide.tex doc/manual.xml
	cd parser.src;make dynare_m.exe
	cd doc;make

$(MZIP_FILE):
	cd ..;zip $(FNAME_MAT) $(DYNARE_ROOT)/matlab/*.m $(DYNARE_ROOT)/matlab/*.exe $(DYNARE_ROOT)/matlab/*.dll $(DYNARE_ROOT)/examples/*.mod $(DYNARE_ROOT)/examples/fs2000/* $(DYNARE_ROOT)/doc/guide.pdf $(DYNARE_ROOT)/doc/manual/*

$(PZIP_FILE): parser.src/dynare_m.exe
	cd ..;zip $(FNAME_PAR) $(DYNARE_ROOT)/parser.src/*.h $(DYNARE_ROOT)/parser.src/*.hh $(DYNARE_ROOT)/parser.src/*.l $(DYNARE_ROOT)/parser.src/*.y $(DYNARE_ROOT)/parser.src/*.c $(DYNARE_ROOT)/parser.src/*.cc $(DYNARE_ROOT)/parser.src/makefile $(DYNARE_ROOT)/parser.src/d.output

install: $(MZIP_FILE) $(PZIP_FILE) 
	cd zips; chmod a+rx $(MZIP_FILE); scp $(MZIP_FILE) pythie.cepremap.cnrs.fr:public_html/mambo/download;chmod a+rx $(PZIP_FILE); scp $(PZIP_FILE) pythie.cepremap.cnrs.fr:public_html/mambo/download

