all:
	(cd rqhmm/src; make clean)
	(cd src; rm -f *.o)
	R CMD INSTALL rqhmm
#	(cd rqhmm_extra/src; make clean)
#	R CMD INSTALL rqhmm_extra

debug:
	(cd rqhmm/src; make clean)
	(cd src; rm -f *.o)
	PKG_CXXFLAGS="-Wall -g -O0" R CMD INSTALL rqhmm
#	(cd rqhmm_extra/src; make clean)
#	PKG_CXXFLAGS="-Wall -g -O0" R CMD INSTALL rqhmm_extra
