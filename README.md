# Louvain-like Methods for Community Detection in Multi-Layer Networks

## Proposed Methods
-EA: Louvain Expansion Method Average (in particular: in EA_s communities are indexed by size and in EA_r are indexed randomly).\
-EVM: Louvain Expansion Method Function F-.\
-EVP: Louvain Expansion Method Function F+.\
-MA: Louvain Multiobjective Method Average (in particular: in MA\_s communities are indexed by size and in MA\_r are indexed randomly).\
-MVM: Louvain Multiobjective Method Function F-.\
-MCP: Louvain Multiobjective Method Function F+.

## State-of-the-art methods (to be added in this folder)
From https://github.com/youweiliang/Multi-view\_Clustering

## Artificial networks 
-adjacent\_matrix\_generator: creates a single-layer graph using the Stochastic Block Model.
-adjacent\_matrix\_generator\_multi: creates a multi-layer graph for the informative case: each layer is informative. Each layer is created by the adjacent\_matrix\_generator.\
-adjacent\_matrix\_generator\_multi\_r: creates a multi-layer graph for the noisy case: first layer noisy and all the other layers informative. Each layer is created by the adjacent_matrix_generator.

## Real Datasets
From https://github.com/melopeo/PM\_SSL/tree/master/realworld\_datasets\
P. Mercado, F. Tudisco, and M. Hein, Generalized Matrix Means for Semi-Supervised Learning with Multilayer Graphs. 

## Evaluation of the final partitions
-confusion_matrix: calculates the confusion matrix. It is used in "wrong" function and input of "NMI" function.\
-reindex_com: reindexes communities. It is used in "confusion_matrix" function.\
-wrong: counts the number of nodes in the wrong community. It is used to calculate the Accuracy of a partition.\
-NMI: calculates the Normalized Mutual Information (NMI) of a partition.

## Tests on Artificial Networks
-run: tests each methods on artificial networks for the informative case.\
-run_n: tests each methods on artificial networks for the noisy case.

## Tests on Real Datasets
-REAL: tests each methods on real datasets.\
-REAL\_n2: tests each methods on real datasets with noise (all informative layers + one noisy layer).\
-REAL\_n3: tests each methods on real datasets with noise (first layer sum of all the real layers, second layer is noise).

## Reference paper
"Louvain-like Methods for Community Detection in Multi-Layer Networks" by Sara Venturini, Andrea Cristofari, Francesco Rinaldi, Francesco Tudisco.

## Authors
- Sara Venturini (e-mail: sara.venturini@math.unipd.it)
- Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
- Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
- Francesco Tudisco (e-mail: francesco.tudisco@gssi.it)
