# Span-Constrained Truss Community Search on Temporal Graphs

# Repository Structure

# 

# data/

# Preprocessed experimental datasets used in our paper.

# 

# ksecforest/

# Source code of the KSECForest index construction and query algorithms.

# 

# ksebgraph/

# Source code of the KSEBGraph index construction and query algorithms.

# 

# Compilation

# 

# Compile the binary executables for both indexes:

# 

# bash build.sh

# 

# 

# After compilation, two executables will be generated:

# 

# ksebgraph

# 

# ksecforest

# 

# Index Construction

# KSEBGraph Index

# ./ksebgraph build <dataset\_name> 0   # TwoPhase algorithm

# ./ksebgraph build <dataset\_name> 1   # OnePass algorithm

# 

# KSECForest Index

# ./ksecforest build <dataset\_name> 0  # TwoPhase algorithm

# ./ksecforest build <dataset\_name> 1  # OnePass algorithm

# 

# 

# After execution, two binary files will be generated:

# 

# A binary file of the original temporal graph

# 

# A binary file of the corresponding index

# 

# Query Processing

# 

# To perform a span-constrained truss community query:

# 

# ./<binary\_name> query <binary\_graph\_file> <binary\_index\_file> <query\_file> <k> <delta>

# 

# 

# Parameters:

# 

# <binary\_name>: ksebgraph or ksecforest

# 

# <binary\_graph\_file>: binary file of the temporal graph

# 

# <binary\_index\_file>: binary file of the constructed index

# 

# <query\_file>: file containing query vertices

# 

# <k>: truss parameter

# 

# <delta>: span constraint

# 

# Notes

# 

# All datasets are preprocessed to match the input format described in the paper.

# 

# Vertices are indexed with consecutive integer IDs.

# 

# The choice between TwoPhase and OnePass controls the index construction strategy.

