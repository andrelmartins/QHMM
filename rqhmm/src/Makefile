LIBSRC:=../../src
SRC:=rqhmm.cpp $(wildcard $(LIBSRC)/*.cpp)
HPP=func_entry.hpp $(wildcard $(LIBSRC)/*.hpp)

# build our MAKEFLAGS
# R needs spaces to be escaped ...
#
MFLGS := "CXXFLAGS=-I../../src -DUSE_RMATH -DRQHMM"
# define a space
empty:=
space:= $(empty) $(empty)
MFLGSESC := $(subst $(space),\$(space),$(MFLGS))


all: rqhmm.so

rqhmm.so: $(SRC) $(HPP)
	MAKEFLAGS=${MFLGSESC} R CMD SHLIB -o $@ $(filter-out %/hmm.cpp, $(SRC))

clean:
	rm -f rqhmm.so $(patsubst %.cpp,%.o,$(SRC))
