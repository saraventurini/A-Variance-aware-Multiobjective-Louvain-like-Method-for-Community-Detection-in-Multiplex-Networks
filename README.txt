This repository contains the codes of the paper  "Louvain-like Methods for Community Detection in Multi-Layer Networks" by Sara Venturini, Andrea Cristofari, Francesco Rinaldi, Francesco Tudisco.

CODES

%Proposed Methods
-EA: Louvain Expansion Method Average 
     In particular: in EA_s communities are indexed by size and in EA_r are indexed randomly
-EVM: Louvain Expansion Method Function F-
-EVP: Louvain Expansion Method Function F+
-MA: Louvain Multiobjective Method Average
     In particular: in MA_s communities are indexed by size and in MA_r are indexed randomly 
-MVM: Louvain Multiobjective Method Function F-
-MCP: Louvain Multiobjective Method Function F+ 

%State-of-the-art methods (to be added in this folder)
From https://github.com/youweiliang/Multi-view_Clustering

%Artificial networks 
-adjacent_matrix_generator: creates a single-layer graph using the Stochastic Block Model 
-adjacent_matrix_generator_multi: creates a multi-layer graph for the informative case: each layer is informative. 
			          Each layer is created by the adjacent_matrix_generator.
-adjacent_matrix_generator_multi_r: creates a multi-layer graph for the noisy case: first layer noisy and all the other layers informative.  
				    Each layer is created by the adjacent_matrix_generator.

%Real Datasets
From https://github.com/melopeo/PM_SSL/tree/master/realworld_datasets
P. Mercado, F. Tudisco, and M. Hein, Generalized Matrix Means for Semi-Supervised Learning with Multilayer Graphs. In NeurIPS 2019.

%Evaluation of the final partitions
-confusion_matrix: calculates the confusion matrix.
		   It is used in "wrong" function and input of "NMI" function.
-reindex_com: reindexes communities. 
	      It is used in "confusion_matrix" function.
-wrong: counts the number of nodes in the wrong community.
	It is used to calculate the Accuracy of a partition.
-NMI: calculates the Normalized Mutual Information (NMI) of a partition 

%Tests on Artificial Networks
-run: tests each methods on artificial networks for the informative case
-run_n: tests each methods on artificial networks for the noisy case

%Tests on Real Datasets
-REAL: tests each methods on real datasets
-REAL_n2: tests each methods on real datasets with noise (all informative layers + one noisy layer)
-REAL_n3: tests each methods on real datasets with noise (first layer sum of all the real layers, second layer is noise)
