all: guide.pdf manual.pdf manual/index.html

guide.pdf: guide.tex
	pdflatex guide

manual/index.html: manual.xml
	xsltproc -stringparam base.dir ./manual/ dynare_html.xsl manual.xml 


manual.pdf: manual.xml
	xsltproc dynare_fo.xsl manual.xml >manual.fo
	pdfxmltex manual.fo	
	pdfxmltex manual.fo	
