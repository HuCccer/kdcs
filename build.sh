#!/bin/bash
set -e
echo "Compiling..."
rm -f kseb ksec

g++ -g -O2 -std=c++17 \
    ./ksebgraph/test/ksebgraphTestMain.cpp \
    ./ksebgraph/src/ksebgraph_mba.cpp \
    -o ./kseb

g++ -g -O2 -std=c++17 \
    ./ksecforest/test/ksecforestTestMain.cpp \
    ./ksecforest/src/ksecforest_mba.cpp \
    -o ./ksec

echo "Compilation finished."
echo
