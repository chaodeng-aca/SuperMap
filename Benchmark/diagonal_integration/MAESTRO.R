library(Seurat)
library(MAESTRO)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

#input data
inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks


#RNA modality
pbmc_rna <- CreateSeuratObject(counts = rna_counts,assay = 'RNA')
pbmc_rna = NormalizeData(object = pbmc_rna)
pbmc_rna = FindVariableFeatures(pbmc_rna)
pbmc_rna = ScaleData(pbmc_rna)
var_gene = VariableFeatures(pbmc_rna)

#ATAC modality
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
pbmc.atac = CreateChromatinAssay(counts = atac_counts,sep = c(":", "-"),annotation = annotations,
                                 fragments = 'pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz')
pbmc.atac = CreateSeuratObject(counts = pbmc.atac,assay = 'ATAC')
pbmc.atac <- RunTFIDF(pbmc.atac)
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = 'q0')
pbmc.atac <- RunSVD(pbmc.atac)


rownames(atac_counts) = chartr(old = '-',new = '_',x = rownames(atac_counts))
MAESTRO_imputation = MAESTRO::ATACCalculateGenescore(inputMat = atac_counts) #calculating gene activity score


rna_seurat = pbmc_rna
atac_seurat = pbmc.atac
ATAC_expression = MAESTRO_imputation ##using the gene activity calculated by maestro

atac_seurat[["ACTIVITY"]] = CreateAssayObject(counts = ATAC_expression)
DefaultAssay(atac_seurat) <- "ACTIVITY"
atac_seurat <- NormalizeData(atac_seurat)
atac_seurat <- ScaleData(atac_seurat, features = rownames(atac_seurat))

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rna_seurat, query = atac_seurat, features =  var_gene,
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")


refdata = (rna_seurat@assays$RNA$data)[var_gene,]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,weight.reduction = atac_seurat[["lsi"]],dims = 2:50)

#transfer label
cell_label = cell_label[colnames(rna_seurat),]
predicted_label <- TransferData(anchorset = transfer.anchors, refdata = cell_label,weight.reduction = atac_seurat[["lsi"]],dims = 2:50)
predicted_label = predicted_label['predicted.id']

#Co-embedding
rna_data = rna_seurat@assays$RNA$data[var_gene,]
atac_data = imputation@data[rownames(rna_data),]
colnames(atac_data) = paste0(colnames(atac_data),'_atac')
coembed_data <- cbind(rna_data,atac_data)
coembed = CreateSeuratObject(coembed_data)
coembed = NormalizeData(coembed)
coembed@assays$RNA$data = coembed_data
coembed <- ScaleData(coembed, features = var_gene, do.scale = FALSE)
coembed <- RunPCA(coembed, features = var_gene, verbose = FALSE)
coembed <- RunUMAP(coembed,dims = 1:50)


