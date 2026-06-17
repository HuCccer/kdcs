# Efficient $(k,\delta)$-Truss Community Search on Dynamic Temporal Graphs

This repository implements the **KSECForest index** for efficient $(k,\delta)$-truss community search on dynamic temporal graphs.

The system supports index construction, query processing, and dynamic maintenance under edge updates.

---

## Repository Structure

- `data/`  
  temporal graph datasets used in experiments.

- `ksecforest/`  
  Core implementation of KSECForest:
  - index construction
  - query processing
  - dynamic maintenance

- `infrastructure/`  
  Basic graph data structures and utilities.

- `ksecforest/test/`  
  Executable programs for construction, query, and update evaluation.

---

## Compilation

Compile all executables using:

```bash
cd kdcs
make -j$(nproc)

After compilation, the following executables will be generated:

constructionTest.exe (binary file) 
addEdgeTest.exe (binary file) 
queryTest.exe (binary file) 


---

## Index Construction
ksecforest/test/constructionTest.exe <graph_file>

After execution, three binary files will be generated:

1. A binary file of the original temporal graph  
2. A binary file of the kspan index
3. A binary file of the ksecforest inndex

---

## Query Processing

To perform a $(k,\delta)$-Truss Community query by k-span index or ksecforest index, run:

ksecforest/test/queryTest.exe <binary temporal_graph_file> <binary kspan_index> <binary ksecforest_index> <query_file> <k_value> <delta_value>

Parameters:

<binary temporal_graph_file> 
Binary file of the temporal graph

<binary kspan_index>
Binary file of the constructed kspan index

<binary ksecforest_index>
Binary file of the constructed ksecforest index

<query_file>  
File containing query vertices

<k_value> 
Truss parameter

<delta_value>
Span constraint

---

## Index Maintenance

To perform ksecforest index maintenance, run:

ksecforest/test/addEdgeTest.exe <binary temporal_graph_file> <binary kspan_index> <binary ksecforest_index> <update_file> 

Parameters:

<binary temporal_graph_file> 
Binary file of the temporal graph

<binary kspan_index>
Binary file of the constructed kspan index

<binary ksecforest_index>
Binary file of the constructed ksecforest index

<update_file>  
list of edges to be inserted


## Example Workflow

# compile
make -j

# construct index
ksecforest/test/constructionTest.exe ../data/askubuntu.txt

# query
ksecforest/test/queryTest.exe ../index/graph_askubuntu.txt ../index/kspan_askubuntu.txt ../index/forest_askubuntu.txt ../query/askubuntu.txt 4 100

# update
ksecforest/test/addEdgeTest.exe ../index/graph_askubuntu.txt ../index/kspan_askubuntu.txt ../index/forest_askubuntu.txt ../update/update_askubuntu.txt 


