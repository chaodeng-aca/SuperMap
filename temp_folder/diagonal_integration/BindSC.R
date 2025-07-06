library(bindSC)
library(umap)
library(Matrix)

inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

#-----data proprecessing-----------------------


#RNA modality
pbmc_rna <- CreateSeuratObject(counts = rna_counts,assay = 'RNA')
pbmc_rna = NormalizeData(object = pbmc_rna)
pbmc_rna = FindVariableFeatures(pbmc_rna)
pbmc_rna = ScaleData(pbmc_rna)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
pbmc.atac = CreateChromatinAssay(counts = atac_counts,sep = c(":", "-"),annotation = annotations,
                                 fragments = 'pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz')
pbmc.atac = CreateSeuratObject(counts = pbmc.atac,assay = 'ATAC')


gas = GeneActivity(object = pbmc.atac,features = var_gene) #calculating gene activity score


atac_dim = dimReduce(dt1 = atac_counts, K=30)
rna = rna_seurat@assays$RNA$data
atac_expression = gas


rna_cluster = cell_label$seurat_clusters
names(rna_cluster) = rownames(cell_label)
atac_cluster = rep(1,ncol(atac_expression))
names(atac_cluster) = colnames(atac_expression)


##run bindsc
bindsc_result <- BiCCA( X = t(rna),
              Y = t(atac_dim), 
              Z0 = t(atac_expression), 
              X.clst = rna_cluster,
              Y.clst = atac_cluster,
              alpha = 0.5, 
              lambda = 0.5,
              K = 10,
              temp.path  = "out",
              num.iteration = 100,
              tolerance = 0.01,
              save = TRUE,
              parameter.optimize = FALSE,
              block.size = 0)


#extraction results
embedding = rbind(bindsc_result$u, bindsc_result$r)
umap_plt <- umap(rbind(bindsc_result$u, bindsc_result$r))
umap_cord = umap_plt$layout


