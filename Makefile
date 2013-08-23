all:
	(cd rqhmm/src; make clean)
	(cd src; rm -f *.o)
	R CMD INSTALL rqhmm

debug:
	(cd rqhmm/src; make clean)
	(cd src; rm -f *.o)
	PKG_CXXFLAGS="-Wall -g -O0" R CMD INSTALL rqhmm
