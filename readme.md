
# _netSmoothLite_: A lite version of _netSmooth_, with only the basic functions

## _netSmooth_ is a network smoothing based method for single cell RNA-seq data imputation 

_netSmooth_ is an R package for network smoothing of single cell RNA sequencing data. Using gene interaction networks such as protein-
protein interactions as priors for gene co-expression, _netsmooth_ improves cell type identification from noisy, sparse scRNA-seq data.
The smoothing method is suitable for other gene-based omics data sets such as proteomics, copy-number variation, etc.

The algorithm uses a network-diffusion based approach which takes in a network (such as PPI network) and gene-expression matrix. The gene 
expression values in the matrix are smoothed using the interaction information in the network. The network-smoothing parameter is 
optimized using a robust clustering approach.

For a detailed exposition, check out [the published paper on F1000Research](https://f1000research.com/articles/7-8/v2).


### Installation of this lite version

	devtools::install_github("sqjin/netSmooth")

### Quick run

	data.impute <- netSmoothLite(data.use, mouse.ppi)
	
	
### How to cite
Please cite the _netSmooth_ paper:

> Ronen J and Akalin A. _netSmooth_: Network-smoothing based imputation for single cell RNA-seq [version 2; referees: 2 approved]. F1000Research 2018, 7:8 (doi: 10.12688/f1000research.13511.2)


