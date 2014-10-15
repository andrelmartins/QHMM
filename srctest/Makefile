all: tests

SOURCES=$(filter-out %hmm.cpp,$(wildcard ../src/*.cpp) $(wildcard *.cpp))

tests: catch.hpp $(SOURCES)
	clang++ -o $@ -Wall -I../src ${SOURCES} -lm

catch.hpp:
	curl https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp > catch.hpp

