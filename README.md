# A Variance-aware Multiobjective Louvain-like Method for Community Detection in Multiplex Networks

## Proposed Methods
- GL: Generalized Louvain (in particular: in GL\_s communities are indexed by size and in GL\_r are indexed randomly).
- EVM: Louvain Expansion Method Function F-.
- EVP: Louvain Expansion Method Function F+.
- MA: Louvain Multiobjective Method Average (in particular: in MA\_s communities are indexed by size and in MA\_r are indexed randomly). 
- MVM: Louvain Multiobjective Method Function F-.
- MVP: Louvain Multiobjective Method Function F+. 

## State-of-the-art methods (to be added in this folder)
### Download
- From https://github.com/youweiliang/Multi-view_Clustering: Co-Regularized Spectral Clustering (CoReg), multi-view clustering via Adaptively Weighted Procrustes (AWP), Multi-view Consensus Graph Clustering (MCGC), and The Power Mean Laplacian Multi-layer Method (PM).
- From https://github.com/cdebacco/MultiTensor: Multitensor expectation maximization method (MT).
- From https://github.com/gilles-didier/MolTi: Subspace Analysis on Grassmann Manifolds (SCML), and Principal Modularity Maximization (PMM).
- From https://github.com/mapequation/infomap: Infomap Information-theoretic generalization of mapequation (IM).
### Help files
- INFOMAP_adjacency_matrix: from adjecency matrix to input file IM method.
- MULTITENSOR_adjacency_matrix: from adjecency matrix to input file for MT method.

## Artificial networks 
### SBM Stochastic Block Model
- adjacent\_matrix\_generator: creates a single-layer graph using the Stochastic Block Model 
- adjacent\_matrix\_generator_multi: creates a multi-layer graph for the informative case: each layer is informative. Each layer is created by the adjacent_matrix_generator.
-adjacent\_matrix\_generator\_multi\_r: creates a multi-layer graph for the noisy case: SOME layer noisy and SOME layers informative. Each layer is created by the adjacent_matrix_generator.
### LFR Lancichinetti-Fortunato-Radicchi
From https://www.santofortunato.net/resources Package 1 \
A. Lancichinetti, S. Fortunato, and F. Radicchi, Benchmark graphs for testing community detection algorithms, Physical review E, vol. 78, no. 4, p. 046110, 2008.

## Real Datasets
From https://github.com/melopeo/PM_SSL/tree/master/realworld_datasets \
P. Mercado, F. Tudisco, and M. Hein, Generalized Matrix Means for Semi-Supervised Learning with Multilayer Graphs. In NeurIPS 2019.

## Evaluation of the final partitions
- confusion_matrix: calculates the confusion matrix. It is used in "wrong" function and input of "NMI" function.
- reindex_com: reindexes communities. It is used in "confusion_matrix" function.
- wrong: counts the number of nodes in the wrong community. It is used to calculate the Accuracy of a partition.
- NMI: calculates the Normalized Mutual Information (NMI) of a partition.

## Tests on Artificial Networks
### SBM
test\_ART: run all the artificial tests SBM (informative and real) reported in the paper.
- run: tests each methods on artificial networks SBM for the informative case.
- run\_n: tests each methods on artificial networks SBM for the noisy case.
### LFR
test_ART_LFR: run all the artificial tests LFR (informative and real) reported in the paper.
- run\_LFR: tests each methods on artificial networks LFR for the informative case .
- run\_n\_LFR: tests each methods on artificial networks LFR for the noisy case.

## Tests on Real Datasets
test\_REAL: run all the tests on real-world networks (informative and real) reported in the paper.
- REAL: tests each methods on real datasets.
- REAL\_n2: tests each methods on real datasets with noise (all informative layers + one noisy layer).
- REAL\_n3: tests each methods on real datasets with noise (first layer sum of all the real layers, second layer is noise).

## Reference paper
"Louvain-like Methods for Community Detection in Multiplex Networks" by Sara Venturini, Andrea Cristofari, Francesco Rinaldi, Francesco Tudisco.

## Authors
- Sara Venturini (e-mail: sara.venturini@math.unipd.it)
- Andrea Cristofari (e-mail: andrea.cristofari@unipd.it)
- Francesco Rinaldi (e-mail: rinaldi@math.unipd.it)
- Francesco Tudisco (e-mail: francesco.tudisco@gssi.it)
