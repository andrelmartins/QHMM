LIBS:= -lgsl -lgslcblas -lm
LIBSRC:=../../src
SRC:=rqhmm.cpp $(wildcard $(LIBSRC)/*.cpp)
HPP=func_entry.hpp $(wildcard $(LIBSRC)/*.hpp)

all: rqhmm.so

rqhmm.so: $(SRC) $(HPP)
	MAKEFLAGS="CXXFLAGS=-I../../src" R CMD SHLIB -o $@ $(filter-out %/hmm.cpp, $(SRC)) $(LIBS)

clean:
	rm -f rqhmm.so $(patsubst %.cpp,%.o,$(SRC))
