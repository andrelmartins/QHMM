

all: logsum_bench fwd_bench

logsum_bench: logsum_bench.cpp ../src/logsum.hpp ../src/logsum.cpp
	$(CXX) -O2 -Wall -o $@ logsum_bench.cpp ../src/logsum.cpp

fwd_bench: fwd_bench.cpp $(wildcard ../src/*.cpp) $(wildcard ../src/*.hpp) $(wildcard ../src/emissions/*.hpp)  $(wildcard ../src/transitions/*.hpp)
	$(CXX) -O2 -g -Wall -o $@ fwd_bench.cpp $(filter-out ../src/hmm.cpp, $(wildcard ../src/*.cpp))
