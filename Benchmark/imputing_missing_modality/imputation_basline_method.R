

#----------------------Signac--------------------------------
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

#input data
inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
pbmc.atac = CreateChromatinAssay(counts = atac_counts,sep = c(":", "-"),annotation = annotations,
                                 fragments = 'pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz')
pbmc.atac = CreateSeuratObject(counts = pbmc.atac,assay = 'ATAC')

gas = GeneActivity(object = pbmc.atac,features = var_gene) #calculating gene activity score


#----------------------MAESTRO--------------------------------

library(MAESTRO)
library(Seurat)

inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks


rownames(atac_counts) = chartr(old = '-',new = '_',x = rownames(atac_counts))
MAESTRO_imputation = MAESTRO::ATACCalculateGenescore(inputMat = atac_counts)


