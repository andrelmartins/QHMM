#!/bin/bash
for m in {1..20}; do x=`echo 16*2^$m | bc`; ./em_test $x | grep -E "^[0123456789]"; done > em_test.txt
