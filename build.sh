#!/bin/bash
set -e
echo "Compiling..."
rm -f ksebgraph ksecforest

g++ -g -O2 -std=c++17 \
    ./ksebgraph/test/ksebgraphTestMain.cpp \
    ./ksebgraph/src/ksebgraph_mba.cpp \
    -o ./ksebgraph

g++ -g -O2 -std=c++17 \
    ./ksecforest/test/ksecforestTestMain.cpp \
    ./ksecforest/src/ksecforest_mba.cpp \
    -o ./ksecforest

echo "Compilation finished."
echo
