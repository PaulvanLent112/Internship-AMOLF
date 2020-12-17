neurons<- readRDS("data/neuron_cds.rds")



#This part is stolen from the CeNGEN BioRxiv_prepint_code_figure3_4A 
L4.neuron <- neurons
# From raw data, this has had all cells with > 10% of UMIs from mitochondrial genes removed.

# The following steps have already been run on the stored monocle object, but I included them here:

# Detecting genes in each cell type
L4.neuron <- detectGenes(L4.neuron)

# Removing genes that are not expressed in any cells
L4.neuron <- L4.neuron[fData(L4.neuron)$num_cells_expressed > 0,]

# This code was used to speed up some of the processing steps (the estimateDispersion and 
# plot_pc_variance_explained steps were particularly time-consuming without it)

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)

# Estimating Size Factors and Dispersion metrics
L4.neuron <- estimateSizeFactors(L4.neuron)
L4.neuron <- estimateDispersions(L4.neuron)

# Selecting the subset of genes to use for PCA
unsup_clustering_genes_L4 <- subset(dispersionTable(L4.neuron), mean_expression > 0.002 & dispersion_empirical > 5)

gene_IDs<- unsup_clustering_genes_L4$gene_id



#Subset genes in their approach for PCA
neurons<-neurons[which(rownames(neurons)%in% gene_IDs),]

dim(neurons)
#Create monocle 3 CDS
exprs<- exprs(neurons)
sample_data<- pData(neurons)
feature_data<- fData(neurons)

library(monocle3)
cds2 <- monocle3::new_cell_data_set(exprs,
                                   cell_metadata= sample_data, gene_metadata=feature_data)


cds2<- preprocess_cds(cds2, method="PCA",num_dim=125, norm_method="log")
cds2<- monocle3::reduce_dimension(cds2, reduction_method="UMAP", max_components = 2, umap.min_dist=0.2, umap.n_neighbors = 75)
cds2<- monocle3::cluster_cells(cds2,resolution=1e-3, k=18)


pData(cds2)$partition<- partitions(cds2)

p1<-plot_cells(cds2)
p2<-plot_cells(cds2, color_cells_by = "Neuron.type", label_groups_by_cluster = FALSE)
p3<-plot_cells(cds2, color_cells_by = "partition", label_groups_by_cluster = FALSE)
p3
