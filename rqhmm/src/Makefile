LIBSRC:=../../src
SRC:=rqhmm.cpp $(wildcard $(LIBSRC)/*.cpp)
HPP=func_entry.hpp $(wildcard $(LIBSRC)/*.hpp)

all: rqhmm.so

rqhmm.so: $(SRC) $(HPP)
	MAKEFLAGS="CXXFLAGS=-I../inst/include" R CMD SHLIB -o $@ $(filter-out %/hmm.cpp, $(SRC))

clean:
	rm -f rqhmm.so $(patsubst %.cpp,%.o,$(SRC))
