all:
	(cd rqhmm/src; make clean)
	(cd src; rm -f *.o)
	R CMD INSTALL rqhmm
