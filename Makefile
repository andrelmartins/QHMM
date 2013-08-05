all:
	(cd rqhmm/src; make clean)
	R CMD INSTALL rqhmm
