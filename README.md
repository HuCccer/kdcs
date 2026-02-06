# Span-Constrained Truss Community Search on Temporal Graphs

This repository contains the implementation and experimental datasets for our study on Span-Constrained Truss Community Search on Temporal Graphs.

---

data: preprocessed experimental datasets used in our paper.

ksecforest: source code of the KSECForest index construction and query algorithms.

ksebgraph: source code of the KSEBGraph index construction and query algorithms.

---

## Compilation

Compile the binary executables for both indexes:

bash build.sh

After compilation, the following executables will be generated:

kseb (binary file) 
ksec (binary file) 

---

## Index Construction

KSEBGraph Index

./kseb build <dataset_name> 0   # TwoPhase algorithm  
./kseb build <dataset_name> 1   # OnePass algorithm  

KSECForest Index

./ksec build <dataset_name> 0  # TwoPhase algorithm  
./ksec build <dataset_name> 1  # OnePass algorithm  

After execution, two binary files will be generated:

1. A binary file of the original temporal graph  
2. A binary file of the corresponding index  

---

## Query Processing

To perform a span-constrained truss community query, run:

./<binary_name> query <binary temporal_graph_file> <binary_index_file> <query_file> <k_value> <delta_value>

Parameters:

<binary_name>  
kseb or ksec

<binary_graph_file>  
Binary file of the temporal graph

<binary_index_file>  
Binary file of the constructed index

<query_file>  
File containing query vertices

<k_value> 
Truss parameter

<delta_value>
 Span constraint

