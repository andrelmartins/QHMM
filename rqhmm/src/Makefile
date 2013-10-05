LIBSRC:=../../src
SRC:=rqhmm.cpp $(wildcard $(LIBSRC)/*.cpp)
HPP=func_entry.hpp $(wildcard $(LIBSRC)/*.hpp)

# define compile options
CXX:=$(shell R CMD config CXX)

OPENMP_FLAG:=
OPENMP_LIB:=

ifeq ($(CXX), g++)
  OPENMP_FLAG := -fopenmp
  OPENMP_LIB := -lgomp
endif

# build our MAKEFLAGS
# R needs spaces to be escaped ...
#
MFLGS := "CXXFLAGS=-I../../src -DUSE_RMATH -DRQHMM $(OPENMP_FLAG)"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))

all: rqhmm.so

rqhmm.so: $(SRC) $(HPP)
	MAKEFLAGS=${MFLGSESC} R CMD SHLIB $(OPENMP_LIB) -o $@ $(filter-out %/hmm.cpp, $(SRC))

clean:
	rm -f rqhmm.so $(patsubst %.cpp,%.o,$(SRC))
