# Span-Constrained Truss Community Search on Temporal Graphs

This repository contains the implementation and experimental datasets for our study on Span-Constrained Truss Community Search on Temporal Graphs.

---

## Repository Structure

.
├── data
└── Preprocessed experimental datasets used in our paper
├── ksecforest
└── Source code of the KSECForest index construction and query algorithms
├── ksebgraph
└── Source code of the KSEBGraph index construction and query algorithms
├── build.sh
└── README.md

**data**

Preprocessed experimental datasets used in our paper.

""ksecforest""

Source code of the KSECForest index construction and query algorithms.

""ksebgraph""

Source code of the KSEBGraph index construction and query algorithms.

---

## Compilation

Compile the binary executables for both indexes:

bash build.sh

After compilation, the following executables will be generated:

ksebgraph  
ksecforest

---

## Index Construction

KSEBGraph Index

./ksebgraph build <dataset_name> 0   # TwoPhase algorithm  
./ksebgraph build <dataset_name> 1   # OnePass algorithm  

KSECForest Index

./ksecforest build <dataset_name> 0  # TwoPhase algorithm  
./ksecforest build <dataset_name> 1  # OnePass algorithm  

After execution, two binary files will be generated:

1. A binary file of the original temporal graph  
2. A binary file of the corresponding index  

---

## Query Processing

To perform a span-constrained truss community query, run:

./<binary_name> query <binary_graph_file> <binary_index_file> <query_file> <k> <delta>

Parameters:

<binary_name>  
ksebgraph or ksecforest

<binary_graph_file>  
Binary file of the temporal graph

<binary_index_file>  
Binary file of the constructed index

<query_file>  
File containing query vertices

<k>  
Truss parameter

<delta>  
Span constraint

