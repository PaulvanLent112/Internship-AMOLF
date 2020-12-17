# Internship-AMOLF
This repository contains:
1. Scripts containing the data analysis for the thesis
2. Additional plots of the gene module expression in CeNGEN
3. The reduced CeNGEN dataset (198 gene modules *52412 cells) 
4. Gene modules
5. Percolation based gene clustering R version (unfinished)



1. Data analysis scripts
Neuron_clustering_createCDS:      From the CeNGEN dataset (20842*52412) to the  average gene module reduced dataset (198*52412)
Neuron_clustering_Downstream1:    UMAP, Differential Expression, Cluster purity, etc.
Reproducing CeNGEN (monocle 3):   PCA based approach as described in (Taylor et. al 2019)
Reference_stm_std:                RGGs (Number of references =300) for different D and N + script
functions                         Some additional functions that were created and are required for some of the data analysis

2. Additional plots
Two directories containing the expression profiles of neuron types or gene modules

3. Reduced CeNGEN dataset
neuron_cds.rds        
The CeNGEN dataset can be downloaded from cengen.org

4. Gene modules
All the gene modules that were found by the algorihtm in .txt files

5. Percolation-based gene clustering R version
- computeclusternumber, computeclustersize, computenoisesizes_maxelements, computenoisesizes5, computingclusterflag, extractallclusters, trackclustersizes.

Contains a few functions that have been transferred from MATLAB to R. The extract_all_cluster is unfinished (in functions folder).


