#!/bin/bash

clang++ -g --coverage -Wall -o fwd_bench fwd_bench.cpp ../src/iter.cpp ../src/logsum.cpp

gcov fwd_bench.cpp | grep -A1 "../../src" | grep -v "\-\-"
