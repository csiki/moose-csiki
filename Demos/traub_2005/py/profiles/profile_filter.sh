#!/bin/bash

for f in *.prof
do
    pprof --text ../gperftools_wrapped.so $f | grep HSolve > $f.log
done
