rm(list = ls())
set.seed(24)
setwd("~")
setwd("/lustre/home/acct-medlqian/medlqian-loop3/data/10X_sc_RNA_Seq/trim66_ko_het_homo_10X_20220527/")

#### load packages ######
.libPaths(new = "/lustre/home/acct-medlqian/medlqian-loop3/R/x86_64-redhat-linux-gnu-library/4.1.3")

# Bioanalysis
suppressPackageStartupMessages(library(Seurat)) # Seurat_3.0.2
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(DoubletFinder)) # https://github.com/chris-mcginnis-ucsf/DoubletFinder
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(monocle3))
# ### install slingshot
# BiocManager::install("slingshot")
suppressPackageStartupMessages(library(slingshot)) # slingshot mannual: https://bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))

suppressPackageStartupMessages(library(cowplot)) # ggplot的辅助插件
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(gridExtra)) # grid.arrange() merge multiple plot in one pdf
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap)) # to add legend in Circos plot
suppressPackageStartupMessages(library(hutils)) # to use Switch, which is the vectorized form of base::switch. like ifelse is vevtorized form of if else. Check: https://www.delftstack.com/howto/r/use-a-vectorized-if-function-with-multiple-conditions-in-r/
suppressPackageStartupMessages(library(grDevices)) 
suppressPackageStartupMessages(library(openxlsx))

#### Step 1: check the data and choose the suitable percent.mt threshold ---------------------------------------
# Load dataset
trim66_het_data <- Read10X(data.dir = "raw_data_merge/MoE_TR1M66-Ko-het_13")
trim66_homo_data <- Read10X(data.dir = "raw_data_merge/MoE_TR1M66-Ko-homo_15")

# filter: remove cells with fewer than 700 or more than 15, 000 UMI counts 
# (criteria used by Michael Greenburg 2017 Nature Neuroscience)
trim66_homo_data_UMI <- trim66_homo_data[, Matrix::colSums(trim66_homo_data)>= 700 & Matrix::colSums(trim66_homo_data) <= 15000]
trim66_het_data_UMI <- trim66_het_data[, Matrix::colSums(trim66_het_data)>= 700 & Matrix::colSums(trim66_het_data) <= 15000]

trim66_het_object <- CreateSeuratObject(counts = trim66_het_data_UMI, 
                                        min.cells = 3, 
                                        min.features = 200, 
                                        project = "het")
trim66_het_object$genotype <- "het"
trim66_het_object

trim66_homo_object <- CreateSeuratObject(counts = trim66_homo_data_UMI, 
                                         min.cells = 3, 
                                         min.features = 200, 
                                         project = "homo") # trim66_homo_object@meta.data$orig.ident
trim66_homo_object$genotype <- "homo"
trim66_homo_object 

trim66_het_object[["percent.mt"]] <- PercentageFeatureSet(object = trim66_het_object, pattern = '^mt-')
trim66_homo_object[["percent.mt"]] <- PercentageFeatureSet(object = trim66_homo_object, pattern = '^mt-')

# nfeature_nCount_percentMt vlnPlot 
pdf("plots_het_homo_raw_data/VlnPlot_trim66_ko_nfeature_nCount_percentMt.pdf", width = 14, height = 7)
VlnPlot(trim66_het_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(trim66_homo_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# nfeature_nCount_percentMt FeaturePlot
pdf("plots_het_homo_raw_data/FeaturePlot_trim66_ko_nfeature_nCount_percentMt.pdf", width = 14, height = 7)
plot1 <- FeatureScatter(trim66_het_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(trim66_het_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
plot1 <- FeatureScatter(trim66_homo_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(trim66_homo_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

##### filter the data (percent.mt < 10)---------------------------------------
# percent.mt < 10
trim66_het_object <- subset(x = trim66_het_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
trim66_homo_object <- subset(x = trim66_homo_object, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

median_nFeature_het <- median(trim66_het_object@meta.data[["nFeature_RNA"]])
median_nFeature_homo <- median(trim66_homo_object@meta.data[["nFeature_RNA"]])

median_nCount_het <- median(trim66_het_object@meta.data[["nCount_RNA"]])
median_nCount_homo <- median(trim66_homo_object@meta.data[["nCount_RNA"]])

save(trim66_homo_object, file = "results_percent.mt_10/trim66_homo_object_filter.rdata")
save(trim66_het_object, file = "results_percent.mt_10/trim66_het_object_filter.rdata")

load("results_percent.mt_10/trim66_homo_object_filter.rdata")
load("results_percent.mt_10/trim66_het_object_filter.rdata")

dim(trim66_homo_object)
# [1] 21981 16662
dim(trim66_het_object)
# [1] 20769 10103

summary(trim66_het_object@meta.data[["nFeature_RNA"]])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201    1987    3064    2742    3628    4890 
summary(trim66_homo_object@meta.data[["nFeature_RNA"]])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201    2632    3366    3037    3796    4913 

summary(trim66_het_object@meta.data[["nCount_RNA"]])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 700    4856    8592    7997   11098   14998
summary(trim66_homo_object@meta.data[["nCount_RNA"]])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 700    7120    9881    9066   11824   15000 

sum(trim66_het_object@meta.data[["nCount_RNA"]])
# 80791274
sum(trim66_homo_object@meta.data[["nCount_RNA"]])
# 151051699

#### Step 2: merge-normalize-findVariable-scale-PCA-harmony #########################################################################################################################################
# harmony to merge the homo and het data
trim66_harmony <- merge(trim66_het_object, y = c(trim66_homo_object))
trim66_harmony <- NormalizeData(trim66_harmony)%>% 
  FindVariableFeatures(mean.function = ExpMean, dispersion.function = LogVMR) %>%
  ScaleData(vars.to.regress = c('percent.mt')) %>% 
  RunPCA(verbose = FALSE)

system.time({trim66_harmony <- RunHarmony(trim66_harmony, group.by.vars = "orig.ident", plot_convergence = TRUE)}) 
# system.time () function will measure how long it takes to run something in R. 

harmony_embeddings <- Embeddings(trim66_harmony, 'harmony')
# In downstream analyses, use the Harmony embeddings instead of PCA.The Harmony algorithm iteratively corrects PCA embeddings.

trim66_harmony <- RunUMAP(trim66_harmony, reduction = "harmony", dims = 1:20)

trim66_harmony <- FindNeighbors(trim66_harmony, reduction = "harmony", dims = 1:20)

save(trim66_harmony, file = "results_percent.mt_10_integrate/trim66_harmony_percent.mt_10_step2.rdata")

#### step 3: choose the best resolution ##############################################################################
load("results_percent.mt_10_integrate/trim66_harmony_percent.mt_10_step2.rdata")


markers_ref <- c("Krt5", "Krt14", "Trp63", "Cxcl14", "Sox2", "Meg3", # HBC
                 "Ascl1", "Sox2", "Lgr5", "Hes6", "Cxcr4", "Sox2", "Ezh2", "Neurog1", "Neurod1", "Tmprss4", "Kit", "Mki67", "Top2a", # GBC
                 "Neurog1", "Neurod1", "Neurod2", "Top2a", "Mki67", "Lhx2", # INP early\mid
                 "Lhx2", "Ebf1", "Gap43", # INP late
                 "Ablim1", "Drd4", "Dbn1", "Dpysl5", "Crmp1", "Ppp2cb", "Marcksl1", "Atf5", "Gnas", "Hdac2", "Gng8", "Dpysl3", "Gap43", "Stmn1", "Stmn2", "Ebf2", "Lhx2", "Cbx8", "Trib3", # immature
                 "Cngb1", "Adcy3", "Slc17a6", "Cnga4", "Cnga2", "Gnal", "Omp", "Gng13", "Stoml3", "Ebf2", "Cbx8", "Rtp1", # mature OSN markers
                 "Cyp2g1", "Cyp1a2", "Cyp2g1", "Notch2", "Notch3", "Hes1", "Hes5", "Gpx6", "Ermn", "Sox2", "Hey1", "Sult1c1", # sus
                 "S100b", "Plp1", "Mpz", "Alx3", # olfactory ensheathing glia
                 "Coch", "Ascl3", "Cftr", "Hepacam2", # olfactory microvillar cells(mv1)
                 "Sox9", "Trpm5", # mv2
                 "Ms4a4c", "Ms4a6c", "Ms4a7", "Pde2a", "Car2", "Cnga3", # Ms4a-expressing chemosensory receptor cell (Ms4)
                 "Foxj1", "Cfap126", "Stoml3", # respiratory ciliated cells 
                 "Sox9", "Scgb1a1", # respiratory gland progenitor cells 
                 "Cyp4b1", "Muc5", "Tff3", # respiratory secretory cells
                 "Tagln", "Myh11", # vascular smooth muscle cell
                 "Sox17", "Eng", # pericyte
                 "Sox9", "Sox10", "Muc5", "Gpx3", # bowman's gland
                 "Cd3d", "Cd3e", "Cd8a", # cd8+ t cell
                 "Cd3d", "Cd3e", "Cd4", "Il7r", # cd4+ t cell
                 "Fgfbp2", "Fcgr3a", "Cx3cr1", # natural killer cell
                 "Cd19", "Cd79a", "Ms4a4a", # b cell
                 "Mzb1", "Sdc1", "Cd79a", # plasma cell
                 "Cd14", "S100a12", "Clec10a", # monocyte
                 "C1qa", "C1qb", "C1qc", # macrophage
                 "Ptprc", "Hba-a1" #blood cell
)#125
markers_ref <- markers_ref[!duplicated(markers_ref)] # 107
markers_ref <- intersect(markers_ref, rownames(trim66_harmony)) # 102

list_FindAllMarkers_top50 <- list()
for (n in 0.1*(3:15)) {
  trim66_harmony <- FindClusters(trim66_harmony, resolution = n) 
  
  trim66_harmony.FindAllMarkers <- FindAllMarkers(object = trim66_harmony, 
                                                  only.pos = TRUE, 
                                                  min.pct = 0.25, 
                                                  logfc.threshold = 0.25)
  
  trim66_harmony.FindAllMarkers_top50 <- trim66_harmony.FindAllMarkers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
  
  list_FindAllMarkers_top50[[(10*n-2)]] <- trim66_harmony.FindAllMarkers_top50
  
  DotPlot_names <- paste0("resolution_", n, ".", "_Dotplot_markers.pdf")
  pdf(width = 125, height = 50, file = DotPlot_names)
  p <- DotPlot(trim66_harmony, features = markers_ref)
  print(p)
  dev.off()
  
  VlnPlot_names = paste0("VlnPlot_Markers_resolution_", n, ".pdf")
  pdf(width = 50, height = 375, file = VlnPlot_names)
  plot <- VlnPlot(trim66_harmony, features = markers_ref, 
                  pt.size = 0, combine = FALSE)
  plot <- CombinePlots(plots = plot, ncol = 1)
  print(plot)
  dev.off()
  
  tsne_names = paste0("tsne_resolution_", n, "_harmony", ".pdf")
  pdf(width = 15, height = 7, file = tsne_names)
  p1 <- DimPlot(trim66_harmony, reduction = "tsne", label = TRUE, repel = TRUE)
  p2 <- DimPlot(trim66_harmony, reduction = "tsne", group.by = "genotype")
  p3 <- p1 + p2
  print(p3)
  dev.off()
  
  umap_names = paste0("umap_resolution_", n, "_harmony", ".pdf")
  pdf(width = 15, height = 7, file = umap_names)
  p4 <- DimPlot(trim66_harmony, reduction = "umap", label = TRUE, repel = TRUE)
  p5 <- DimPlot(trim66_harmony, reduction = "umap", group.by = "genotype")
  p6 <- p4 + p5
  print(p6)
  dev.off()
  
} 

n = length(list_FindAllMarkers_top50)
n
names(list_FindAllMarkers_top50) = paste0("resolution_", 0.1*(3:15))
write.xlsx(list_FindAllMarkers_top50, file = "./FindAllMarkers.resolution_0.3_1.5_percent.mt_10.xlsx")

# list_featurePlot_tsne <- list()
# pdf("featurePlot_markers_ref_tsne.pdf", width = 7, height = 7)
# for (i in c(1:102))
# {list_featurePlot_tsne[[i]] <- FeaturePlot(object = trim66_harmony, 
#                                            features = markers_ref[i], 
#                                            reduction = "tsne")
# }
# print(list_featurePlot_tsne)
# dev.off()

list_featurePlot_umap <- list()
pdf("featurePlot_markers_ref_umap.pdf", width = 7, height = 7)
for (i in c(1:102))
{list_featurePlot_umap[[i]] <- FeaturePlot(object = trim66_harmony, 
                                           features = markers_ref[i], 
                                           reduction = "umap")
}
print(list_featurePlot_umap)
dev.off()

#### step 4: add metadata (resolution = 0.9) --------------------------------------------------------
load("results_percent.mt_10_integrate/trim66_harmony_percent.mt_10_step2.rdata")

trim66_harmony <- FindClusters(trim66_harmony, resolution = 0.9) # 0-24 cluster

save(trim66_harmony, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_step3.rdata")

### create idents and add to metadata
# celltype
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_step3.rdata")
Idents(trim66_harmony) <- "seurat_clusters"
trim66_harmony <- RenameIdents(trim66_harmony, 
                               '0' = 'mOSN', '1' = 'mOSN', '2' = 'mOSN', '3' = 'mOSN', '5' = 'mOSN', '6' = 'mOSN', '18' = 'mOSN', # 7
                               '4' = 'imOSN', '7' = 'imOSN', '11' = 'imOSN', # 3
                               '12' = 'INP', 
                               '22' = 'GBC', 
                               '21' = 'HBC', 
                               '8' = 'SUS', '13' = 'SUS', 
                               '15' = 'MV', '17' = 'MV', # microvilliar
                               '10' = 'Ms4a', '14' = 'Ms4a', '20' = 'Ms4a', 
                               '19' = 'T cell', 
                               '23' = 'B cell', 
                               '9' = 'Blood cell', '16' = 'Blood cell', 
                               '24' = 'UD')
trim66_harmony$celltype <- Idents(trim66_harmony)
Idents(trim66_harmony) <- "celltype"
trim66_harmony@active.ident <- factor(trim66_harmony@active.ident, levels = c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD"))

# celltype.genotype
trim66_harmony$celltype.genotype <- paste(Idents(trim66_harmony), trim66_harmony$genotype, sep = "_")

# celltype_number
Idents(trim66_harmony) <- "seurat_clusters"
trim66_harmony <- RenameIdents(trim66_harmony, 
                               '0' = 'mOSN_1', '1' = 'mOSN_2', '2' = 'mOSN_3', '3' = 'mOSN_4', '5' = 'mOSN_5', '6' = 'mOSN_6', '18' = 'mOSN_7', 
                               '4' = 'imOSN_1', '7' = 'imOSN_2', '11' = 'imOSN_3', 
                               '12' = 'INP', 
                               '22' = 'GBC', 
                               '21' = 'HBC', 
                               '8' = 'SUS_1', '13' = 'SUS_2', 
                               '15' = 'MV_1', '17' = 'MV_2', # microvilliar
                               '10' = 'Ms4a_1', '14' = 'Ms4a_2', '20' = 'Ms4a_3', 
                               '19' = 'T cell', 
                               '23' = 'B cell', 
                               '9' = 'Blood cell_1', '16' = 'Blood cell_2', 
                               '24' = 'UD')
trim66_harmony$celltype_number <- Idents(trim66_harmony)
trim66_harmony@active.ident <- factor(trim66_harmony@active.ident, levels = c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                                                                              "imOSN_1", "imOSN_2", "imOSN_3", 
                                                                              "INP", 
                                                                              "GBC", 
                                                                              "HBC", 
                                                                              "SUS_1", "SUS_2", 
                                                                              "MV_1", "MV_2", 
                                                                              "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                                                                              "T cell", 
                                                                              "B cell", 
                                                                              "Blood cell_1", "Blood cell_2", 
                                                                              "UD"))
trim66_harmony$celltype_number.genotype <- paste(Idents(trim66_harmony), trim66_harmony$genotype, sep = "_")

# seurat_clusters
Idents(trim66_harmony) <- "seurat_clusters"
trim66_harmony[["seurat_clusters"]]$seurat_clusters <- factor(trim66_harmony[["seurat_clusters"]]$seurat_clusters, levels = c(0, 1, 2, 3, 5, 6, 18, 4, 7, 11, 12, 22, 21, 8, 13, 15, 17, 10, 14, 20, 19, 23, 9, 16, 24))

save(trim66_harmony, file = 'trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')

load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')

colnames(trim66_harmony@meta.data)
# metadata
# [1] "orig.ident" "nCount_RNA" "nFeature_RNA" 
# [4] "genotype" "percent.mt" "RNA_snn_res.0.9" 
# [7] "seurat_clusters" "celltype" "celltype.genotype" 
# [10] "celltype_number" "celltype_number.genotype"

##### plot 4.1: tsne, umap, vlnplot, dotplot, marker_gene_featurePlot ####
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
markers_plot_markers <-c("Omp", "Gnal", "Cnga2", # mature OSN markers
                         "Gap43", "Gng8", # immature
                         "Neurog1", "Neurod1", # INP
                         "Mki67", # GBC
                         "Krt5", "Krt14", # HBC
                         "Cyp2g1", "Cyp1a2", # SUS
                         "Coch", "Hepacam2", "Sox9", # mv
                         "Ms4a4c", "Ms4a6c", "Ms4a7", # Ms4a
                         "Cd3d", # T cell
                         "Cd79a", # B cell
                         "Ptprc", "Hba-a1" #blood cell
) # 22个 

### umap
pdf(width = 8, height = 7, file = 'results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/umapPlot_with_celltype_percent.mt_10_resolution_0.9.pdf')
p <- DimPlot(trim66_harmony, reduction = "umap", group.by = "celltype", label = TRUE, 
             cols = c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2", "cyan3", "palegreen3", "peachpuff3", "lightslateblue", "plum3", "palevioletred3", "sienna3", "azure4"), 
             pt.size = 1)
print(p)
dev.off()

### umap split by genotype
Idents(trim66_harmony) <- "genotype"
pdf(width = 13, height = 7, file = 'results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/umapPlot_splitByGenotype_percent.mt_10_resolution_0.9.pdf')
p <- DimPlot(trim66_harmony, reduction = "umap", split.by = "genotype", label = TRUE, 
             cols = c("springgreen3", "lightcoral"), pt.size = 1)
print(p)
dev.off()

### umap group by genotype
# pdf(width = 8, height = 7, file = 'results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/umapPlot_groupByGenotype_percent.mt_10_resolution_0.9.pdf')
# list <- list()
# for (i in 1:9) {
#   
#   list[[i]] <- DimPlot(trim66_harmony, reduction = "umap", group.by = "genotype", label = TRUE, 
#                        cols = c("springgreen3", alpha("lightcoral", i/10)), 
#                        pt.size = 1)
#   
# }
# print(list)
# dev.off()
# # choose the third plot alpha = 0.3

### umap labeled by celltype and labeled by genotype 
pdf(width = 15, height = 7, file = 'results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/umapPlot_percent.mt_10_resolution_0.9.pdf')
p1 <- DimPlot(trim66_harmony, reduction = "umap", group.by = "celltype", label = TRUE, 
              cols = c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2", "cyan3", "palegreen3", "peachpuff3", "lightslateblue", "plum3", "palevioletred3", "sienna3", "azure4"), 
              pt.size = 1)
p2 <- DimPlot(trim66_harmony, reduction = "umap", group.by = "genotype", label = TRUE, 
              cols = c("springgreen3", alpha("lightcoral", 0.3)), 
              pt.size = 1)
p <- p1+p2
print(p)
dev.off()

### violin plot without x_axis_text according to the order of celltype 
pdf(file = "results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/VlnPlot_selected_markers_percent.mt_10_resolution_0.9.pdf", width = 25, height = 40)
plot <- VlnPlot(trim66_harmony, features = markers_plot_markers, 
                group.by = "seurat_clusters", 
                cols = c(rep("royalblue3", 7), rep("cornflowerblue", 3), rep("lightblue3", 1), rep("paleturquoise2", 1), rep("cyan3", 1), 
                         rep("palegreen3", 2), rep("peachpuff3", 2), rep("lightslateblue", 3), rep("plum3", 1), rep("palevioletred3", 1), 
                         rep("palevioletred3", 1), rep("sienna3", 2), rep("azure4", 1)), 
                pt.size = 0, ncol = 1, combine = T)
plot <- plot & 
  theme(axis.title.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank()) & 
  xlab(NULL) & 
  ylab(NULL)
print(plot)
dev.off()

### violin plot with cluster name
pdf(file = "results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/VlnPlot_selected_markers_percent.mt_10_resolution_0.9_with_xlab.pdf", width = 25, height = 40)
plot <- VlnPlot(trim66_harmony, features = markers_plot_markers, 
                group.by = "seurat_clusters", 
                cols = c(rep("royalblue3", 7), rep("cornflowerblue", 3), rep("lightblue3", 1), rep("paleturquoise2", 1), rep("cyan3", 1), 
                         rep("palegreen3", 2), rep("peachpuff3", 2), rep("lightslateblue", 3), rep("plum3", 1), rep("palevioletred3", 1), 
                         rep("palevioletred3", 1), rep("sienna3", 2), rep("azure4", 1)), 
                pt.size = 0, ncol = 1, combine = T)
plot <- plot & xlab(NULL) & ylab(NULL)
print(plot)
dev.off()

### violin plot with celltype_number name
Idents(trim66_harmony) <- "celltype_number"
pdf(width = 25, height = 60, file = "results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/VlnPlot_with_celltype_number_selected_markers_percent.mt_10_resolution_0.9.pdf")
plot <- VlnPlot(trim66_harmony, features = markers_plot_markers, 
                group.by = "celltype_number", 
                cols = c(rep("royalblue3", 7), rep("cornflowerblue", 3), rep("lightblue3", 1), rep("paleturquoise2", 1), rep("cyan3", 1), 
                         rep("palegreen3", 2), rep("peachpuff3", 2), rep("lightslateblue", 3), rep("plum3", 1), rep("palevioletred3", 1), 
                         rep("palevioletred3", 1), rep("sienna3", 2), rep("azure4", 1)), 
                pt.size = 0, combine = T, ncol = 1)
print(plot)
dev.off()

### dotplot with cluster
pdf(width = 12, height = 10, file = "results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/dotplot_with_cluster_selected_markers_percent.mt_10_resolution_0.9.pdf")
p <- DotPlot(trim66_harmony, features = rev(markers_plot_markers), group.by = "seurat_clusters", cols = c("lightgrey", "dodgerblue3")) + 
  RotatedAxis()+
  coord_flip () 
print(p)
dev.off()

### dotplot with celltype_number name
pdf(width = 12, height = 10, file = "results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/dotplot_with_celltype_selected_markers_percent.mt_10_resolution_0.9.pdf")
p <- DotPlot(trim66_harmony, features = rev(markers_plot_markers), group.by = "celltype_number", cols = c("lightgrey", "dodgerblue3")) + 
  RotatedAxis()+
  coord_flip () 
print(p)
dev.off()

### featurePlot_umap
list_featurePlot_umap <- list()
pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/featurePlot_markers_ref_umap.pdf", width = 7, height = 7)
for (i in c(1:22))
{list_featurePlot_umap[[i]] <- FeaturePlot(object = trim66_harmony, 
                                           features = markers_plot_markers[i], 
                                           reduction = "umap")
}
print(list_featurePlot_umap)
dev.off()
##### step 4.2.1: 查看 celltype/cluster 的细胞数和比例 ---------------------------------
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
# celltype
celltype <- c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
celltype_data <- list()
celltype_subset <- list()
for (i in 1:length(celltype)) {
  Idents(trim66_harmony) <- "celltype"
  celltype_subset[[i]] <- subset(x = trim66_harmony, idents = celltype[i])
  celltype_data[[i]] <- t(as.data.frame(dim(celltype_subset[[i]])))
}

multirbind <-function(dat = list(), ...){
  if(length(dat) < 2)return(as.data.frame(dat))
  rbinddat <-dat[[1]]
  dat[[1]] <-NULL
  for(i in dat){
    rbinddat <-rbind(rbinddat, i)
  }
  return(rbinddat)
}

celltype_pattern <- multirbind(celltype_data)
rownames(celltype_pattern) <- celltype 
celltype_pattern <- as.data.frame(celltype_pattern)
colnames(celltype_pattern) <- c("gene", "cellNumber")
celltype_pattern$percent <- celltype_pattern[, "cellNumber"]/ncol(trim66_harmony)
write.csv(celltype_pattern, "celltype_data_pattern.csv")

#celltype_number
celltype_number <- c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                     "imOSN_1", "imOSN_2", "imOSN_3", 
                     "INP", 
                     "GBC", 
                     "HBC", 
                     "SUS_1", "SUS_2", 
                     "MV_1", "MV_2", 
                     "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                     "T cell", 
                     "B cell", 
                     "Blood cell_1", "Blood cell_2", 
                     "UD")
celltype_number_data <- list()
celltype_number_subset <- list()
for (i in 1:25) {
  Idents(trim66_harmony) <- "celltype_number"
  celltype_number_subset[[i]] <- subset(x = trim66_harmony, idents = celltype_number[i])
  celltype_number_data[[i]] <- t(as.data.frame(dim(celltype_number_subset[[i]])))
}
celltype_number_pattern <- multirbind(celltype_number_data)
rownames(celltype_number_pattern) <- celltype_number
colSums(celltype_number_pattern)
celltype_number_pattern <- as.data.frame(celltype_number_pattern)
colnames(celltype_number_pattern) <- c("gene", "cellNumber")
celltype_number_pattern$percent <- celltype_number_pattern[, "cellNumber"]/ncol(trim66_harmony)
write.csv(celltype_number_pattern, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/celltype_number_data_pattern.csv")

##### step 4.2.2：查看 het/homo celltype/celltype_nunmber 细胞数和比例 ------------------------------
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "genotype"
het_trim66_data <- subset(x = trim66_harmony, idents = "het") # 10103
homo_trim66_data <- subset(x = trim66_harmony, idents = "homo") # 16662

# het_celltype
het_celltype <- c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
het_celltype_data <- list()
het_celltype_subset <- list()
for (i in 1:12) {
  Idents(het_trim66_data) <- "celltype"
  het_celltype_subset[[i]] <- subset(x = het_trim66_data, idents = het_celltype[i])
  het_celltype_data[[i]] <- t(as.data.frame(dim(het_celltype_subset[[i]])))
}

het_celltype_pattern <- multirbind(het_celltype_data)
rownames(het_celltype_pattern) <- het_celltype
colSums(het_celltype_pattern)
het_celltype_pattern <- as.data.frame(het_celltype_pattern)
colnames(het_celltype_pattern) <- c("gene", "cellNumber")
het_celltype_pattern$percent <- het_celltype_pattern[, "cellNumber"]/ncol(het_trim66_data)
write.csv(het_celltype_pattern, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_celltype_data_pattern.csv")

# het_celltype_number
het_celltype_number <- c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                         "imOSN_1", "imOSN_2", "imOSN_3", 
                         "INP", 
                         "GBC", 
                         "HBC", 
                         "SUS_1", "SUS_2", 
                         "MV_1", "MV_2", 
                         "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                         "T cell", 
                         "B cell", 
                         "Blood cell_1", "Blood cell_2", 
                         "UD")
het_celltype_number_data <- list()
het_celltype_number_subset <- list()
for (i in 1:25) {
  Idents(het_trim66_data) <- "celltype_number"
  het_celltype_number_subset[[i]] <- subset(x = het_trim66_data, idents = het_celltype_number[i])
  het_celltype_number_data[[i]] <- t(as.data.frame(dim(het_celltype_number_subset[[i]])))
}

het_celltype_number_pattern <- multirbind(het_celltype_number_data)
rownames(het_celltype_number_pattern) <- het_celltype_number
colSums(het_celltype_number_pattern)
het_celltype_number_pattern <- as.data.frame(het_celltype_number_pattern)
colnames(het_celltype_number_pattern) <- c("gene", "cellNumber")
het_celltype_number_pattern$percent <- het_celltype_number_pattern[, "cellNumber"]/ncol(het_trim66_data)
write.csv(het_celltype_number_pattern, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_celltype_number_data_pattern.csv")

# homo_celltype
homo_celltype <- c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
homo_celltype_data <- list()
homo_celltype_subset <- list()
for (i in 1:12) {
  Idents(homo_trim66_data) <- "celltype"
  homo_celltype_subset[[i]] <- subset(x = homo_trim66_data, idents = homo_celltype[i])
  homo_celltype_data[[i]] <- t(as.data.frame(dim(homo_celltype_subset[[i]])))
}


homo_celltype_pattern <- multirbind(homo_celltype_data)
rownames(homo_celltype_pattern) <- homo_celltype
colSums(homo_celltype_pattern)
homo_celltype_pattern <- as.data.frame(homo_celltype_pattern)
colnames(homo_celltype_pattern) <- c("gene", "cellNumber")
homo_celltype_pattern$percent <- homo_celltype_pattern[, "cellNumber"]/ncol(homo_trim66_data)
write.csv(homo_celltype_pattern, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_celltype_data_pattern.csv")

# homo_celltype_number
homo_celltype_number <- c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                          "imOSN_1", "imOSN_2", "imOSN_3", 
                          "INP", 
                          "GBC", 
                          "HBC", 
                          "SUS_1", "SUS_2", 
                          "MV_1", "MV_2", 
                          "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                          "T cell", 
                          "B cell", 
                          "Blood cell_1", "Blood cell_2", 
                          "UD")
homo_celltype_number_data <- list()
homo_celltype_number_subset <- list()
for (i in 1:25) {
  Idents(homo_trim66_data) <- "celltype_number"
  homo_celltype_number_subset[[i]] <- subset(x = homo_trim66_data, idents = homo_celltype_number[i])
  homo_celltype_number_data[[i]] <- t(as.data.frame(dim(homo_celltype_number_subset[[i]])))
}


homo_celltype_number_pattern <- multirbind(homo_celltype_number_data)
rownames(homo_celltype_number_pattern) <- homo_celltype_number
colSums(homo_celltype_number_pattern)
homo_celltype_number_pattern <- as.data.frame(homo_celltype_number_pattern)
colnames(homo_celltype_number_pattern) <- c("gene", "cellNumber")
homo_celltype_number_pattern$percent <- homo_celltype_number_pattern[, "cellNumber"]/ncol(homo_trim66_data)
write.csv(homo_celltype_number_pattern, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_celltype_number_data_pattern.csv")

### integrate the data of cell distribution
# celltype
celltype_all <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/celltype_data_pattern.csv", header = T, row.names = 1)
colnames(celltype_all) <- paste0("all_", colnames(celltype_all))

celltype_het <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_celltype_data_pattern.csv", header = T, row.names = 1)
colnames(celltype_het) <- paste0("het_", colnames(celltype_het))

celltype_homo <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_celltype_data_pattern.csv", header = T, row.names = 1)
colnames(celltype_homo) <- paste0("homo_", colnames(celltype_homo))

celltype_data_pattern <- cbind(celltype_all, celltype_het) %>% cbind(celltype_homo)
write.csv(celltype_data_pattern, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/summary_celltype_data_pattern.csv")

# celltype_number
celltype_number_all <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/celltype_number_data_pattern.csv", header = T, row.names = 1)
colnames(celltype_number_all) <- paste0("all_", colnames(celltype_number_all))

celltype_number_het <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_celltype_number_data_pattern.csv", header = T, row.names = 1)
colnames(celltype_number_het) <- paste0("het_", colnames(celltype_number_het))

celltype_number_homo <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_celltype_number_data_pattern.csv", header = T, row.names = 1)
colnames(celltype_number_homo) <- paste0("homo_", colnames(celltype_number_homo))

celltype_number_data_pattern <- cbind(celltype_number_all, celltype_number_het) %>% cbind(celltype_number_homo)
write.csv(celltype_number_data_pattern, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/summary_celltype_number_data_pattern.csv")

##### plot 4.2: barplot for het/homo celltype distribution ######
# cell_distribution_pattern_in_het_and_homo.csv merged by excel
cell_distribution <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/cell_distribution_pattern_in_het_and_homo.csv", header = T, row.names = 1, encoding = "UTF-8")
cell_distribution_barplot_data <- cell_distribution %>% 
  rownames_to_column(var = 'celltype') %>% 
  pivot_longer( cols = c("Het_percent", "Homo_percent"), 
                names_to = 'sampletype', 
                values_to = 'percent') %>% 
  mutate(celltype = factor(celltype, level = c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/celltype_distribtion_pattern_in_het_and_homo.pdf", height = 7, width = 7)
ggplot(data = cell_distribution_barplot_data, aes(sampletype, percent)) +
  geom_bar(aes(color = celltype, fill = celltype), stat = "identity", position = "stack") +
  scale_fill_manual(values = c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2", "cyan3", "palegreen3", "peachpuff3", "lightslateblue", "plum3", "palevioletred3", "sienna3", "azure4")) + 
  scale_colour_manual(values = c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2", "cyan3", "palegreen3", "peachpuff3", "lightslateblue", "plum3", "palevioletred3", "sienna3", "azure4")) + 
  theme_bw() + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0))
dev.off()

#### step 5: DE analysis ----------------------------------------------------------------
##### step 5.1: DE analysis (celltype) ####
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
celltype_subset <- list()
avg_celltype <- list()
for (i in 1:12){
  celltype <-c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
  Idents(trim66_harmony) <- "celltype"
  celltype_subset[[i]] <- subset(trim66_harmony, idents = celltype[i])
  Idents(celltype_subset[[i]]) <- "genotype"
  avg_celltype[[i]] <- as.data.frame(log1p(AverageExpression(celltype_subset[[i]], verbose = FALSE)$RNA))
  avg_celltype[[i]]$gene <- rownames(avg_celltype[[i]])
}
save(avg_celltype, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_celltype.rdata")

# DE analysis results (default filtering parameters)
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_celltype.rdata")

DEGs_MAST <- list()

avg_DEGs_MAST <- avg_celltype

for (i in 1:12){
  celltype <-c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
  celltype_het <- paste0(celltype, "_het")
  celltype_homo <- paste0(celltype, "_homo")
  Idents(trim66_harmony) <- "celltype.genotype"
  DEGs_MAST[[i]] <- as.data.frame(FindMarkers(trim66_harmony, ident.1 = celltype_homo[i] , ident.2 = celltype_het[i], test.use = "MAST", verbose = FALSE))
  avg_DEGs_MAST[[i]]$p_adj_MAST <- DEGs_MAST[[i]][match(rownames(avg_DEGs_MAST[[i]]), rownames( DEGs_MAST[[i]])), ]$p_val_adj
  avg_DEGs_MAST[[i]]$avg_logFC_MAST <- DEGs_MAST[[i]][match(rownames(avg_DEGs_MAST[[i]]), rownames(DEGs_MAST[[i]])), ]$avg_log2FC
  avg_DEGs_MAST[[i]]$pct.1 <- DEGs_MAST[[i]][match(rownames(avg_DEGs_MAST[[i]]), rownames(DEGs_MAST[[i]])), ]$pct.1
  avg_DEGs_MAST[[i]]$pct.2 <- DEGs_MAST[[i]][match(rownames(avg_DEGs_MAST[[i]]), rownames(DEGs_MAST[[i]])), ]$pct.2
}

save(DEGs_MAST, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/DEGs_MAST.rdata")
save(avg_DEGs_MAST, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST.rdata")

celltype <-c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
names(avg_DEGs_MAST) = paste0(celltype, "_MAST")
write.xlsx(avg_DEGs_MAST, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/avg_DEGs_MAST.xlsx")
names(DEGs_MAST) = paste0(celltype, "_MAST")
write.xlsx(DEGs_MAST, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/DEGs_MAST.xlsx")

# DE analysis (no filter)
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_celltype.rdata")

DEGs_MAST_all <- list()

avg_DEGs_MAST_all <- avg_celltype

for (i in 1:12){
  celltype <-c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
  celltype_het <- paste0(celltype, "_het")
  celltype_homo <- paste0(celltype, "_homo")
  Idents(trim66_harmony) <- "celltype.genotype"
  DEGs_MAST_all[[i]] <- as.data.frame(FindMarkers(trim66_harmony, ident.1 = celltype_homo[i] , ident.2 = celltype_het[i] , test.use = "MAST", verbose = FALSE, 
                                                  logfc.threshold = 0, min.pct = 0))
  avg_DEGs_MAST_all[[i]]$p_adj_MAST <- DEGs_MAST_all[[i]][match(rownames(avg_DEGs_MAST_all[[i]]), rownames( DEGs_MAST_all[[i]])), ]$p_val_adj
  avg_DEGs_MAST_all[[i]]$avg_logFC_MAST <- DEGs_MAST_all[[i]][match(rownames(avg_DEGs_MAST_all[[i]]), rownames(DEGs_MAST_all[[i]])), ]$avg_log2FC
  avg_DEGs_MAST_all[[i]]$pct.1 <- DEGs_MAST_all[[i]][match(rownames(avg_DEGs_MAST_all[[i]]), rownames(DEGs_MAST_all[[i]])), ]$pct.1
  avg_DEGs_MAST_all[[i]]$pct.2 <- DEGs_MAST_all[[i]][match(rownames(avg_DEGs_MAST_all[[i]]), rownames(DEGs_MAST_all[[i]])), ]$pct.2
}
save(DEGs_MAST_all, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/DEGs_MAST_all.rdata")

save(avg_DEGs_MAST_all, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST_all.rdata")

celltype <- c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD")
names(avg_DEGs_MAST_all) = paste0(celltype, "_MAST")
write.xlsx(avg_DEGs_MAST_all, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST_all.xlsx")
names(DEGs_MAST_all) = paste0(celltype, "_MAST")
write.xlsx(DEGs_MAST_all, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/DEGs_MAST_all.xlsx")

###### Plot 5.1.1: scatter plot (no filter) (celltype) ####
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST.rdata")
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST_all.rdata")

DEG_list_MAST <- list()
for (i in 1:12) {
  DEG_list_MAST[[i]] <- subset(avg_DEGs_MAST[[i]], (p_adj_MAST < 0.05 & abs(avg_logFC_MAST)> 0.585) == TRUE)
}

for(i in c(1:12)){
  print(nrow(DEG_list_MAST[[i]] ))
}
# [1] 42
# [1] 25
# [1] 5
# [1] 4
# [1] 2
# [1] 153
# [1] 10
# [1] 23
# [1] 8
# [1] 89
# [1] 307
# [1] 1

for (i in c(1:12)) {
  scattername <-paste0("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/DE_analysis/", c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD"), "_", "scatterPlot.pdf") 
  pdf(scattername[i], width = 7, height = 7)
  p1 <- ggplot(avg_DEGs_MAST[[i]], aes(x = het, y = homo)) + 
    geom_point(colour = "grey", size = 1) +
    xlab("het")+
    ylab("homo")+
    geom_point(data = subset(avg_DEGs_MAST[[i]], (p_adj_MAST < 0.05 & avg_logFC_MAST> 0.585) == TRUE), colour = "springgreen3") +
    geom_point(data = subset(avg_DEGs_MAST[[i]], (p_adj_MAST < 0.05 & avg_logFC_MAST < -0.585) == TRUE), colour = "steelblue3") +
    ggtitle(" p_adj_MAST < 0.05 and abs(log2FC_MAST)>0.585") + 
    theme(plot.title = element_text(size = 8)) +
    theme_bw()+
    theme(axis.line = element_line(colour = "black"), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank())
  p2 <- LabelPoints(plot = p1, points = rownames(avg_DEGs_MAST[[i]][(avg_DEGs_MAST[[i]]$p_adj_MAST < 0.05) & (abs(avg_DEGs_MAST[[i]]$avg_logFC_MAST)> 0.585), ]), repel = TRUE)
  print(p2)
  dev.off()
}

DEGs_OlfrOrTaar_mOSNOrimOSN <- list()
for (i in 1:2) {
  DEGs_OlfrOrTaar_mOSNOrimOSN[[i]] <- DEG_list_MAST[[i]][grep(pattern = "^Olfr|^Taar", rownames(DEG_list_MAST[[i]])), ]
}
celltype <- c("mOSN", "imOSN")
names(DEGs_OlfrOrTaar_mOSNOrimOSN) = celltype
write.xlsx(DEGs_OlfrOrTaar_mOSNOrimOSN, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/DEGs_OlfrOrTaar_mOSNOrimOSN.xlsx")

###### Plot 5.1.2: volcanol Plot (no filter) (celltype) #####
load("avg_DEGs_MAST_all.rdata")

for (i in c(1:12)) {
  # add a column $diffexpressed
  avg_DEGs_MAST_all[[i]]$diffexpressed <- 0
  avg_DEGs_MAST_all[[i]]$diffexpressed[avg_DEGs_MAST[[i]]$avg_logFC_MAST > 0.585 & avg_DEGs_MAST[[i]]$p_adj_MAST < 0.05] <- "log2FC>0.585 & FDR < 0.05"
  avg_DEGs_MAST_all[[i]]$diffexpressed[avg_DEGs_MAST[[i]]$avg_logFC_MAST < -0.585 & avg_DEGs_MAST[[i]]$p_adj_MAST < 0.05] <- "log2FC <-0.585 & FDR < 0.05"
  # add a column $delabel
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  avg_DEGs_MAST_all[[i]]$delabel <- NA
  avg_DEGs_MAST_all[[i]]$delabel[avg_DEGs_MAST_all[[i]]$diffexpressed != 0] <- avg_DEGs_MAST_all[[i]]$gene[avg_DEGs_MAST_all[[i]]$diffexpressed != 0]
  avg_DEGs_MAST_all[[i]]$p_adj_MAST[which(avg_DEGs_MAST_all[[i]]$p_adj_MAST < 10^(-300))] <- 10^(-300)
  volcanoname <-paste0("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/DE_analysis/", c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a", "T cell", "B cell", "Blood cell", "UD"), "_", "volcanoPlot.pdf") 
  pdf(volcanoname[i], width = 7, height = 7)
  p1 <- ggplot(avg_DEGs_MAST_all[[i]], aes(x = avg_logFC_MAST, y = -log10(p_adj_MAST), col = diffexpressed, label = delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values = c( "black", "steelblue3", "springgreen3")) +
    labs(x = "log2FC", y = "-log10(FDR)")+
    theme_bw()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(title = NULL))+
    theme(legend.position = "none") 
  print(p1)
  dev.off()
}

##### step 5.2: DE analysis (cluster)--------------------------------------------------------
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')

celltype_number_subset <- list()
avg_celltype_number <- list()
for (i in 1:25){
  celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                      "imOSN_1", "imOSN_2", "imOSN_3", 
                      "INP", 
                      "GBC", 
                      "HBC", 
                      "SUS_1", "SUS_2", 
                      "MV_1", "MV_2", 
                      "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                      "T cell", 
                      "B cell", 
                      "Blood cell_1", "Blood cell_2", 
                      "UD")
  Idents(trim66_harmony) <- "celltype_number"
  celltype_number_subset[[i]] <- subset(trim66_harmony, idents = celltype_number[i])
  Idents(celltype_number_subset[[i]]) <- "genotype"
  avg_celltype_number[[i]] <- as.data.frame(log1p(AverageExpression(celltype_number_subset[[i]], verbose = FALSE)$RNA))
  avg_celltype_number[[i]]$gene <- rownames(avg_celltype_number[[i]])
}
save(avg_celltype_number, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_celltype_number.rdata")

# DE analysis (default filtering parameters)
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_celltype_number.rdata")

DEGs_MAST_cluster <- list()

avg_DEGs_MAST_cluster <- avg_celltype_number

for (i in 1:25){
  celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                      "imOSN_1", "imOSN_2", "imOSN_3", 
                      "INP", 
                      "GBC", 
                      "HBC", 
                      "SUS_1", "SUS_2", 
                      "MV_1", "MV_2", 
                      "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                      "T cell", 
                      "B cell", 
                      "Blood cell_1", "Blood cell_2", 
                      "UD")
  celltype_number_het <- paste0(celltype_number, "_het")
  celltype_number_homo <- paste0(celltype_number, "_homo")
  Idents(trim66_harmony) <- "celltype_number.genotype"
  DEGs_MAST_cluster[[i]] <- as.data.frame(FindMarkers(trim66_harmony, ident.1 = celltype_number_homo[i], ident.2 = celltype_number_het[i] , test.use = "MAST", verbose = FALSE))
  avg_DEGs_MAST_cluster[[i]]$p_adj_MAST <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames( DEGs_MAST_cluster[[i]])), ]$p_val_adj
  avg_DEGs_MAST_cluster[[i]]$avg_logFC_MAST <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames(DEGs_MAST_cluster[[i]])), ]$avg_log2FC
  avg_DEGs_MAST_cluster[[i]]$pct.1 <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames(DEGs_MAST_cluster[[i]])), ]$pct.1
  avg_DEGs_MAST_cluster[[i]]$pct.2 <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames(DEGs_MAST_cluster[[i]])), ]$pct.2
}

save(DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/DEGs_MAST_cluster.rdata")
save(avg_DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST_cluster.rdata")

celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                    "imOSN_1", "imOSN_2", "imOSN_3", 
                    "INP", 
                    "GBC", 
                    "HBC", 
                    "SUS_1", "SUS_2", 
                    "MV_1", "MV_2", 
                    "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                    "T cell", 
                    "B cell", 
                    "Blood cell_1", "Blood cell_2", 
                    "UD")
names(avg_DEGs_MAST_cluster) = paste0(celltype_number, "_MAST")
write.xlsx(avg_DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/avg_DEGs_MAST_cluster.xlsx")

celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                    "imOSN_1", "imOSN_2", "imOSN_3", 
                    "INP", 
                    "GBC", 
                    "HBC", 
                    "SUS_1", "SUS_2", 
                    "MV_1", "MV_2", 
                    "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                    "T cell", 
                    "B cell", 
                    "Blood cell_1", "Blood cell_2", 
                    "UD")
names(DEGs_MAST_cluster) = paste0(celltype_number, "_MAST")
write.xlsx(DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/DEGs_MAST_cluster.xlsx")

load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_celltype_number.rdata")

DEGs_MAST_cluster <- list()

avg_DEGs_MAST_cluster <- avg_celltype_number

for (i in 1:25){
  celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                      "imOSN_1", "imOSN_2", "imOSN_3", 
                      "INP", 
                      "GBC", 
                      "HBC", 
                      "SUS_1", "SUS_2", 
                      "MV_1", "MV_2", 
                      "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                      "T cell", 
                      "B cell", 
                      "Blood cell_1", "Blood cell_2", 
                      "UD")
  celltype_number_het <- paste0(celltype_number, "_het")
  celltype_number_homo <- paste0(celltype_number, "_homo")
  Idents(trim66_harmony) <- "celltype_number.genotype"
  DEGs_MAST_cluster[[i]] <- as.data.frame(FindMarkers(trim66_harmony, ident.1 = celltype_number_homo[i], ident.2 = celltype_number_het[i] , test.use = "MAST", verbose = FALSE, 
                                                      logfc.threshold = 0, min.pct = 0))
  avg_DEGs_MAST_cluster[[i]]$p_adj_MAST <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames( DEGs_MAST_cluster[[i]])), ]$p_val_adj
  avg_DEGs_MAST_cluster[[i]]$avg_logFC_MAST <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames(DEGs_MAST_cluster[[i]])), ]$avg_log2FC
  avg_DEGs_MAST_cluster[[i]]$pct.1 <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames(DEGs_MAST_cluster[[i]])), ]$pct.1
  avg_DEGs_MAST_cluster[[i]]$pct.2 <- DEGs_MAST_cluster[[i]][match(rownames(avg_DEGs_MAST_cluster[[i]]), rownames(DEGs_MAST_cluster[[i]])), ]$pct.2
}
save(DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/DEGs_MAST_cluster_all.rdata")
save(avg_DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST_cluster_all.rdata")

celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                    "imOSN_1", "imOSN_2", "imOSN_3", 
                    "INP", 
                    "GBC", 
                    "HBC", 
                    "SUS_1", "SUS_2", 
                    "MV_1", "MV_2", 
                    "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                    "T cell", 
                    "B cell", 
                    "Blood cell_1", "Blood cell_2", 
                    "UD")
names(avg_DEGs_MAST_cluster) = paste0(celltype_number, "_MAST")
write.xlsx(avg_DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/avg_DEGs_MAST_cluster_all.xlsx")


celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                    "imOSN_1", "imOSN_2", "imOSN_3", 
                    "INP", 
                    "GBC", 
                    "HBC", 
                    "SUS_1", "SUS_2", 
                    "MV_1", "MV_2", 
                    "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                    "T cell", 
                    "B cell", 
                    "Blood cell_1", "Blood cell_2", 
                    "UD")
names(DEGs_MAST_cluster) = paste0(celltype_number, "_MAST")
write.xlsx(DEGs_MAST_cluster, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/DEGs_MAST_cluster_all.xlsx")

###### Plot 5.2.1: scatter plot (no filter) (celltype_number) ####
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST_cluster_all.rdata")

DEG_list_MAST <- list()

for (i in 1:25) {
  DEG_list_MAST[[i]] <- subset(avg_DEGs_MAST_cluster[[i]], (p_adj_MAST < 0.05 & abs(avg_logFC_MAST)> 0.585) == TRUE)
}

for(i in c(1:25)){
  print(nrow(DEG_list_MAST[[i]] ))
}
# [1] 54
# [1] 57
# [1] 38
# [1] 48
# [1] 67
# [1] 108
# [1] 35
# [1] 4
# [1] 66
# [1] 5
# [1] 5
# [1] 4
# [1] 2
# [1] 158
# [1] 299
# [1] 4
# [1] 13
# [1] 11
# [1] 7
# [1] 53
# [1] 8
# [1] 89
# [1] 451
# [1] 11
# [1] 1

for (i in c(1:25)) {
  scattername <-paste0("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/DE_analysis/",
                       c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                         "imOSN_1", "imOSN_2", "imOSN_3", 
                         "INP", 
                         "GBC", 
                         "HBC", 
                         "SUS_1", "SUS_2", 
                         "MV_1", "MV_2", 
                         "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                         "T cell", 
                         "B cell", 
                         "Blood cell_1", "Blood cell_2", 
                         "UD"), "_", "scatterPlot.pdf") 
  pdf(scattername[i], width = 7, height = 7)
  p1 <- ggplot(avg_DEGs_MAST_cluster[[i]], aes(x = het, y = homo)) + 
    geom_point(colour = "grey", size = 1) +
    geom_point(data = subset(avg_DEGs_MAST_cluster[[i]], (p_adj_MAST < 0.05 & abs(avg_logFC_MAST)> 0.585) == TRUE), colour = "red") +
    ggtitle(" p_adj_MAST < 0.05 and abs(log2FC_MAST)>0.585") + theme(plot.title = element_text(size = 8)) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  p2 <- LabelPoints(plot = p1, points = rownames(avg_DEGs_MAST_cluster[[i]][(avg_DEGs_MAST_cluster[[i]]$p_adj_MAST < 0.05) & (abs(avg_DEGs_MAST_cluster[[i]]$avg_logFC_MAST)> 0.585), ]), repel = TRUE)
  print(p2)
  dev.off()
}

DEGs_OlfrOrTaar_cluster <- list()
for (i in 1:10) {
  DEGs_OlfrOrTaar_cluster[[i]] <- DEG_list_MAST[[i]][grep(pattern = "^Olfr|^Taar", rownames(DEG_list_MAST[[i]])), ]
}
celltype_number <-c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                    "imOSN_1", "imOSN_2", "imOSN_3")
names(DEGs_OlfrOrTaar_cluster) = celltype_number
write.xlsx(DEGs_OlfrOrTaar_cluster, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/DEGs_OlfrOrTaar_cluster.xlsx")

###### Plot 5.2.2: volcanol Plot (no filter) (celltype_number) #####
load("results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/avg_DEGs_MAST_cluster_all.rdata")
avg_DEGs_MAST_cluster_all <- avg_DEGs_MAST_cluster

for (i in c(1:25)) {
  # add a column $diffexpressed
  avg_DEGs_MAST_cluster_all[[i]]$diffexpressed <- 0
  avg_DEGs_MAST_cluster_all[[i]]$diffexpressed[avg_DEGs_MAST_cluster[[i]]$avg_logFC_MAST > 0.585 & avg_DEGs_MAST_cluster[[i]]$p_adj_MAST < 0.05] <- "log2FC>0.585 & FDR < 0.05"
  avg_DEGs_MAST_cluster_all[[i]]$diffexpressed[avg_DEGs_MAST_cluster[[i]]$avg_logFC_MAST < -0.585 & avg_DEGs_MAST_cluster[[i]]$p_adj_MAST < 0.05] <- "log2FC <-0.585 & FDR < 0.05"
  # add a column $delabel
  # Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
  avg_DEGs_MAST_cluster_all[[i]]$delabel <- NA
  avg_DEGs_MAST_cluster_all[[i]]$delabel[avg_DEGs_MAST_cluster_all[[i]]$diffexpressed != 0] <- avg_DEGs_MAST_cluster_all[[i]]$gene[avg_DEGs_MAST_cluster_all[[i]]$diffexpressed != 0]
  volcanoname <-paste0("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/DE_analysis/",
                       c("mOSN_1", "mOSN_2", "mOSN_3", "mOSN_4", "mOSN_5", "mOSN_6", "mOSN_7", 
                         "imOSN_1", "imOSN_2", "imOSN_3", 
                         "INP", 
                         "GBC", 
                         "HBC", 
                         "SUS_1", "SUS_2", 
                         "MV_1", "MV_2", 
                         "Ms4a_1", "Ms4a_2", "Ms4a_3", 
                         "T cell", 
                         "B cell", 
                         "Blood cell_1", "Blood cell_2", 
                         "UD"), "_", "volcanoPlot.pdf") 
  pdf(volcanoname[i], width = 7, height = 7)
  p1 <- ggplot(avg_DEGs_MAST_cluster_all[[i]], aes(x = avg_logFC_MAST, y = -log10(p_adj_MAST), col = diffexpressed, label = delabel)) + 
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    scale_color_manual(values = c( "black", "lightskyroyalblue3", "seagreen4")) +
    labs(x = "log2FC", y = "-log10(FDR)")+
    theme_bw()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill = guide_legend(title = NULL))+
    theme(legend.position = "none") 
  print(p1)
  dev.off()
}

#### Step 6: UMI threshold for receptor expression (only_functional_receptors) ####
##### Plot 6.1: mOSN: UMI threshold: 1 2 3 4 8 16 32 64 128 256 过滤后 0/single/multiple 占比曲线 (only_functional_receptors)---------
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"
# OR info
all_gene_information <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", sep = "\t", header = T)
OR_and_taar_functional_reference <- all_gene_information$gene_name[intersect(grep(pattern = "protein_coding", all_gene_information$gene_type), grep(pattern = "^Olfr|^Taar", all_gene_information$gene_name))]
OR_and_taar_functional_reference <- OR_and_taar_functional_reference[!duplicated(OR_and_taar_functional_reference)] # 1152

###### het mOSN #####
het_mOSN <- subset(x = trim66_harmony, idents = "mOSN_het")
het_mOSN_functional_Olfr_Taar_expression <- het_mOSN[OR_and_taar_functional_reference, ] # 1036 
dim(het_mOSN_functional_Olfr_Taar_expression)
# UMI matrix
het_mOSN_functional_Olfr_Taar_UMI <- GetAssayData(object = het_mOSN_functional_Olfr_Taar_expression, slot = "counts") #(eg. “counts”, “data”, or “scale.data”)
het_mOSN_functional_Olfr_Taar_UMI <- as.matrix(het_mOSN_functional_Olfr_Taar_UMI)
# > threshold c(0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)
counts_filtered_het_mOSN <- list()
cellNumber_multiple_OR_het_mOSN <- list()
cellNumber_0_OR_het_mOSN <- list()
threshold <- c(0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)
threshold <- as.numeric(threshold)
for (i in threshold) {
  counts_filtered_het_mOSN[[which(threshold == i)]] <- het_mOSN_functional_Olfr_Taar_UMI[, colSums(het_mOSN_functional_Olfr_Taar_UMI > i)>0]
  cellNumber_multiple_OR_het_mOSN[[which(threshold == i)]] <- as.data.frame(table(colSums(counts_filtered_het_mOSN[[which(threshold == i)]] > i)))
  colnames(cellNumber_multiple_OR_het_mOSN[[which(threshold == i)]]) <- c("number", paste0("het_mOSN_", i))
  cellNumber_0_OR_het_mOSN[[which(threshold == i)]] <- ncol(het_mOSN_functional_Olfr_Taar_UMI)-ncol(counts_filtered_het_mOSN[[which(threshold == i)]])
}
cellNumber_0_OR_het_mOSN <- as.numeric(cellNumber_0_OR_het_mOSN)

multifull_join <-function(dat = list(), ...){
  if(length(dat) < 2)return(as.data.frame(dat))
  full_joindat <-dat[[1]]
  dat[[1]] <-NULL
  for(i in dat){
    full_joindat <-full_join(full_joindat, i, by = colnames(i)[1])
  }
  return(full_joindat)
}

cellNumber_multiple_OR_het_mOSN_all_threshold <- multifull_join(cellNumber_multiple_OR_het_mOSN)
cellNumber_multiple_OR_het_mOSN_all_threshold[is.na(cellNumber_multiple_OR_het_mOSN_all_threshold)] <- 0
cellNumber_multiple_OR_het_mOSN_all_threshold_t <- t(cellNumber_multiple_OR_het_mOSN_all_threshold)
cellNumber_multiple_OR_het_mOSN_all_threshold_t <- as.data.frame(apply(cellNumber_multiple_OR_het_mOSN_all_threshold_t, 2, as.numeric))
rownames(cellNumber_multiple_OR_het_mOSN_all_threshold_t) <- colnames(cellNumber_multiple_OR_het_mOSN_all_threshold)
colnames(cellNumber_multiple_OR_het_mOSN_all_threshold_t) <- cellNumber_multiple_OR_het_mOSN_all_threshold_t[1, ]
cellNumber_multiple_OR_het_mOSN_all_threshold_t$multiple <- rowSums(cellNumber_multiple_OR_het_mOSN_all_threshold_t[, 2:11])
cellNumber_multiple_OR_het_mOSN_all_threshold_t$"0" <- c(0, cellNumber_0_OR_het_mOSN)
cellNumber_multiple_OR_het_mOSN_all_threshold_t <- cellNumber_multiple_OR_het_mOSN_all_threshold_t[-1, ]
cellNumber_multiple_OR_het_mOSN_all_threshold_0_1_mutiple <- cellNumber_multiple_OR_het_mOSN_all_threshold_t[, c(1, 12, 13)]
cellNumber_multiple_OR_het_mOSN_all_threshold_0_1_mutiple$threshold <- threshold

cellNumber_multiple_OR_het_mOSN_all_threshold_0_1_mutiple_longer <- cellNumber_multiple_OR_het_mOSN_all_threshold_0_1_mutiple %>% 
  pivot_longer( cols = c("0", "1", "multiple"), 
                names_to = 'OR_number', 
                values_to = 'count') %>% 
  mutate(percent = count/ncol(het_mOSN), 
         threshold = threshold + 1,
         OR_number = factor(OR_number, levels = c("0", "1", "multiple")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/linePlot_het_mOSN_functional_Olfr_Taar_UMI_thresold_percent_change_xAxis.pdf", width = 8, height = 7)
ggplot(data = cellNumber_multiple_OR_het_mOSN_all_threshold_0_1_mutiple_longer, aes(x = as.factor(threshold), y = percent, group = OR_number, color = OR_number)) +
  geom_line()+
  scale_color_manual(values = c("azure4", "deepskyblue2", "hotpink2"))+
  geom_point()+
  xlab("threshold (# of UMIs)")+#横坐标名称
  ylab("Percent of mOSNs")+#纵坐标名称
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())
dev.off()

###### homo mOSN #####
homo_mOSN <- subset(x = trim66_harmony, idents = "mOSN_homo")
homo_mOSN_functional_Olfr_Taar_expression <- homo_mOSN[OR_and_taar_functional_reference, ] # 1036 
# UMI matrix
homo_mOSN_functional_Olfr_Taar_UMI <- GetAssayData(object = homo_mOSN_functional_Olfr_Taar_expression, slot = "counts") #(eg. “counts”, “data”, or “scale.data”). 
homo_mOSN_functional_Olfr_Taar_UMI <- as.matrix(homo_mOSN_functional_Olfr_Taar_UMI)
# > threshold c(0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)
counts_filtered_homo_mOSN <- list()
cellNumber_multiple_OR_homo_mOSN <- list()
cellNumber_0_OR_homo_mOSN <- list()
threshold <- c(0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)
threshold <- as.numeric(threshold)
for (i in threshold) {
  counts_filtered_homo_mOSN[[which(threshold == i)]] <- homo_mOSN_functional_Olfr_Taar_UMI[, colSums(homo_mOSN_functional_Olfr_Taar_UMI > i)>0]
  cellNumber_multiple_OR_homo_mOSN[[which(threshold == i)]] <- as.data.frame(table(colSums(counts_filtered_homo_mOSN[[which(threshold == i)]] > i)))
  colnames(cellNumber_multiple_OR_homo_mOSN[[which(threshold == i)]]) <- c("number", paste0("homo_mOSN_", i))
  cellNumber_0_OR_homo_mOSN[[which(threshold == i)]] <- ncol(homo_mOSN_functional_Olfr_Taar_UMI)-ncol(counts_filtered_homo_mOSN[[which(threshold == i)]])
}
cellNumber_0_OR_homo_mOSN <- as.numeric(cellNumber_0_OR_homo_mOSN)

cellNumber_multiple_OR_homo_mOSN_all_threshold <- multifull_join(cellNumber_multiple_OR_homo_mOSN)
cellNumber_multiple_OR_homo_mOSN_all_threshold[is.na(cellNumber_multiple_OR_homo_mOSN_all_threshold)] = 0
cellNumber_multiple_OR_homo_mOSN_all_threshold_t <- t(cellNumber_multiple_OR_homo_mOSN_all_threshold)
cellNumber_multiple_OR_homo_mOSN_all_threshold_t <- as.data.frame(apply(cellNumber_multiple_OR_homo_mOSN_all_threshold_t, 2, as.numeric))
rownames(cellNumber_multiple_OR_homo_mOSN_all_threshold_t) <- colnames(cellNumber_multiple_OR_homo_mOSN_all_threshold)
colnames(cellNumber_multiple_OR_homo_mOSN_all_threshold_t) <- cellNumber_multiple_OR_homo_mOSN_all_threshold_t[1, ]
cellNumber_multiple_OR_homo_mOSN_all_threshold_t$multiple <- rowSums(cellNumber_multiple_OR_homo_mOSN_all_threshold_t[, 2:11])
cellNumber_multiple_OR_homo_mOSN_all_threshold_t$"0" <- c(0, cellNumber_0_OR_homo_mOSN)
cellNumber_multiple_OR_homo_mOSN_all_threshold_t <- cellNumber_multiple_OR_homo_mOSN_all_threshold_t[-1, ]
cellNumber_multiple_OR_homo_mOSN_all_threshold_0_1_mutiple <- cellNumber_multiple_OR_homo_mOSN_all_threshold_t[, which(colnames(cellNumber_multiple_OR_homo_mOSN_all_threshold_t) %in% c("1", "0", "multiple"))]
cellNumber_multiple_OR_homo_mOSN_all_threshold_0_1_mutiple$threshold <- threshold

cellNumber_multiple_OR_homo_mOSN_all_threshold_0_1_mutiple_longer <- cellNumber_multiple_OR_homo_mOSN_all_threshold_0_1_mutiple %>% 
  pivot_longer( cols = c("0", "1", "multiple"), 
                names_to = 'OR_number', 
                values_to = 'count') %>% 
  mutate(percent = count/ncol(homo_mOSN), 
         threshold = threshold + 1,
         OR_number = factor(OR_number, levels = c("0", "1", "multiple")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/linePlot_homo_mOSN_functional_Olfr_Taar_UMI_thresold_percent_change_xAxis.pdf", width = 8, height = 7)
ggplot(data = cellNumber_multiple_OR_homo_mOSN_all_threshold_0_1_mutiple_longer, aes(x = as.factor(threshold), y = percent, group = OR_number, color = OR_number)) +
  geom_line()+
  scale_color_manual(values = c("azure4", "deepskyblue2", "hotpink2"))+
  geom_point()+
  xlab("threshold (# of UMIs)")+#横坐标名称
  ylab("Percent of mOSNs")+#纵坐标名称
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank()) 
dev.off()

##### Plot 6.2: imOSN: UMI threshold: 1 2 3 4 8 16 32 64 128 256 过滤后 0/single/multiple 占比曲线 (only_functional_receptors)---------
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"
# OR info
all_gene_information <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", sep = "\t", header = T)
OR_and_taar_functional_reference <- all_gene_information$gene_name[intersect(grep(pattern = "protein_coding", all_gene_information$gene_type), grep(pattern = "^Olfr|^Taar", all_gene_information$gene_name))]
OR_and_taar_functional_reference <- OR_and_taar_functional_reference[!duplicated(OR_and_taar_functional_reference)] # 1152

###### het imOSN #####
# UMI matrix
het_imOSN <- subset(x = trim66_harmony, idents = "imOSN_het")
het_imOSN_functional_Olfr_Taar_expression <- het_imOSN[OR_and_taar_functional_reference, ] # 1036 

# UMI > threshold (0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)
het_imOSN_functional_Olfr_Taar_UMI <- as.matrix(GetAssayData(object = het_imOSN_functional_Olfr_Taar_expression, slot = "counts")) # (eg. “counts”, “data”, or “scale.data”)
counts_filtered_het_imOSN <- list()
cellNumber_multiple_OR_het_imOSN <- list()
cellNumber_0_OR_het_imOSN <- list()
threshold <- c(0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)

for (i in threshold) {
  # UMI>i过滤后, 筛选有表达受体的细胞，获得表达矩阵
  counts_filtered_het_imOSN[[which(threshold == i)]] <- as.data.frame(het_imOSN_functional_Olfr_Taar_UMI[, colSums(het_imOSN_functional_Olfr_Taar_UMI > i)>0])
  # UMI>i判断后, colsum计算每个细胞表达的几个受体，table整合
  cellNumber_multiple_OR_het_imOSN[[which(threshold == i)]] <- as.data.frame(table(colSums(counts_filtered_het_imOSN[[which(threshold == i)]] > i)))
  # 修改colname
  colnames(cellNumber_multiple_OR_het_imOSN[[which(threshold == i)]]) <- c("number", paste0("het_imOSN_", i))
  # 增加不表达受体的细胞数
  cellNumber_0_OR_het_imOSN[[which(threshold == i)]] <- ncol(het_imOSN_functional_Olfr_Taar_UMI)-ncol(counts_filtered_het_imOSN[[which(threshold == i)]])
}
cellNumber_0_OR_het_imOSN <- as.numeric(cellNumber_0_OR_het_imOSN)

multifull_join <-function(dat = list(), ...){
  if(length(dat) < 2)return(as.data.frame(dat))
  full_joindat <-dat[[1]]
  dat[[1]] <-NULL
  for(i in dat){
    full_joindat <-full_join(full_joindat, i, by = colnames(i)[1])
  }
  return(full_joindat)
}

cellNumber_multiple_OR_het_imOSN_all_threshold <- multifull_join(cellNumber_multiple_OR_het_imOSN)

cellNumber_multiple_OR_het_imOSN_all_threshold_0_1_multiple <- cellNumber_multiple_OR_het_imOSN_all_threshold %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  dplyr::select(-1) %>% 
  t() %>%
  as.data.frame() %>% 
  mutate("threshold" = threshold, 
         "multiple" = rowSums(dplyr::select(., 2:6)), 
         '0' = cellNumber_0_OR_het_imOSN, 
         '1' = V1) %>%
  dplyr::select('0', '1', 'multiple', 'threshold')


# plot
cellNumber_multiple_OR_het_imOSN_all_threshold_0_1_mutiple_longer <- cellNumber_multiple_OR_het_imOSN_all_threshold_0_1_multiple %>% 
  pivot_longer( cols = c("0", "1", "multiple"), 
                names_to = 'OR_number', 
                values_to = 'count') %>% 
  mutate(percent = count/ncol(het_imOSN), 
         threshold = threshold + 1,
         OR_number = factor(OR_number, levels = c("0", "1", "multiple")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/linePlot_het_imOSN_functional_Olfr_Taar_UMI_thresold_percent_change_xAxis.pdf", width = 8, height = 7)
ggplot(data = cellNumber_multiple_OR_het_imOSN_all_threshold_0_1_mutiple_longer, aes(x = as.factor(threshold), y = percent, group = OR_number, color = OR_number)) +
  geom_line()+
  scale_color_manual(values = c("azure4", "deepskyblue2", "hotpink2"))+
  geom_point()+
  xlab("threshold (# of UMIs)")+#横坐标名称
  ylab("Percent of mOSNs")+#纵坐标名称
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())
dev.off()

###### homo imOSN #####
# UMI matix
homo_imOSN <- subset(x = trim66_harmony, idents = "imOSN_homo")
homo_imOSN_functional_Olfr_Taar_expression <- homo_imOSN[OR_and_taar_functional_reference, ] # 1036 
dim(homo_imOSN_functional_Olfr_Taar_expression)
# UMI > threshold (0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)
homo_imOSN_functional_Olfr_Taar_UMI <- as.matrix(GetAssayData(object = homo_imOSN_functional_Olfr_Taar_expression, slot = "counts")) # (eg. “counts”, “data”, or “scale.data”)
counts_filtered_homo_imOSN <- list()
cellNumber_multiple_OR_homo_imOSN <- list()
cellNumber_0_OR_homo_imOSN <- list()
threshold <- c(0, 1, 2, 3, 4, 7, 15, 31, 63, 127, 255)

for (i in threshold) {
  # UMI>i过滤后, 筛选有表达受体的细胞，获得表达矩阵
  counts_filtered_homo_imOSN[[which(threshold == i)]] <- as.data.frame(homo_imOSN_functional_Olfr_Taar_UMI[, colSums(homo_imOSN_functional_Olfr_Taar_UMI > i)>0])
  # UMI>i判断后, colsum计算每个细胞表达的几个受体，table整合
  cellNumber_multiple_OR_homo_imOSN[[which(threshold == i)]] <- as.data.frame(table(colSums(counts_filtered_homo_imOSN[[which(threshold == i)]] > i)))
  # 修改colname
  colnames(cellNumber_multiple_OR_homo_imOSN[[which(threshold == i)]]) <- c("number", paste0("homo_imOSN_", i))
  # 增加不表达受体的细胞数
  cellNumber_0_OR_homo_imOSN[[which(threshold == i)]] <- ncol(homo_imOSN_functional_Olfr_Taar_UMI)-ncol(counts_filtered_homo_imOSN[[which(threshold == i)]])
}
cellNumber_0_OR_homo_imOSN <- as.numeric(cellNumber_0_OR_homo_imOSN)

multifull_join <-function(dat = list(), ...){
  if(length(dat) < 2)return(as.data.frame(dat))
  full_joindat <-dat[[1]]
  dat[[1]] <-NULL
  for(i in dat){
    full_joindat <-full_join(full_joindat, i, by = colnames(i)[1])
  }
  return(full_joindat)
}

cellNumber_multiple_OR_homo_imOSN_all_threshold <- multifull_join(cellNumber_multiple_OR_homo_imOSN)

cellNumber_multiple_OR_homo_imOSN_all_threshold_0_1_multiple <- cellNumber_multiple_OR_homo_imOSN_all_threshold %>%
  mutate_all(~ifelse(is.na(.), 0, .)) %>%
  dplyr::select(-1) %>% 
  t() %>%
  as.data.frame() %>% 
  mutate("threshold" = threshold, 
         "multiple" = rowSums(dplyr::select(., 2:6)), 
         '0' = cellNumber_0_OR_homo_imOSN, 
         '1' = V1) %>%
  dplyr::select('0', '1', 'multiple', 'threshold')


# plot
cellNumber_multiple_OR_homo_imOSN_all_threshold_0_1_mutiple_longer <- cellNumber_multiple_OR_homo_imOSN_all_threshold_0_1_multiple %>% 
  pivot_longer( cols = c("0", "1", "multiple"), 
                names_to = 'OR_number', 
                values_to = 'count') %>% 
  mutate(percent = count/ncol(homo_imOSN), 
         threshold = threshold + 1,
         OR_number = factor(OR_number, levels = c("0", "1", "multiple")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/linePlot_homo_imOSN_functional_Olfr_Taar_UMI_thresold_percent_change_xAxis.pdf", width = 8, height = 7)
ggplot(data = cellNumber_multiple_OR_homo_imOSN_all_threshold_0_1_mutiple_longer, aes(x = as.factor(threshold), y = percent, group = OR_number, color = OR_number)) +
  geom_line()+
  scale_color_manual(values = c("azure4", "deepskyblue2", "hotpink2"))+
  geom_point()+
  xlab("threshold (# of UMIs)")+#横坐标名称
  ylab("Percent of mOSNs")+#纵坐标名称
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())
dev.off()
#### Step 7: 受体表达矩阵 (only_functional_receptors) ####
##### Step 7.1: mOSN (UMI > 1) het/homo single/multiple 受体表达矩阵 (only_functional_receptors)---------------------------------------------------------
# 1: normalized expression matrix 
# 2: single/mutiple het/homo normalized expression matrix

load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"
# OR info
all_gene_information <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", sep = "\t", header = T)
OR_and_taar_functional_reference <- all_gene_information$gene_name[intersect(grep(pattern = "protein_coding", all_gene_information$gene_type), grep(pattern = "^Olfr|^Taar", all_gene_information$gene_name))]
OR_and_taar_functional_reference <- OR_and_taar_functional_reference[!duplicated(OR_and_taar_functional_reference)] # 1152

### mOSN_het
het_mOSN <- subset(x = trim66_harmony, idents = "mOSN_het")
het_mOSN_functional_Olfr_Taar_expression <- het_mOSN[OR_and_taar_functional_reference, ] # 1036 
# normalized expression matrix
het_mOSN_normalized_expression <- GetAssayData(object = het_mOSN_functional_Olfr_Taar_expression, slot = "data") #(eg. “counts”, “data”, or “scale.data”). 
het_mOSN_normalized_expression <- as.matrix(het_mOSN_normalized_expression)
# UMI matrix
het_mOSN_UMI <- GetAssayData(object = het_mOSN_functional_Olfr_Taar_expression, slot = "counts") #(eg. “counts”, “data”, or “scale.data”). 
het_mOSN_UMI <- as.matrix(het_mOSN_UMI)
all(colnames(het_mOSN_normalized_expression) == colnames(het_mOSN_UMI))
# filter (UMI > 1)
het_mOSN_normalized_expression_filtered <- het_mOSN_normalized_expression
het_mOSN_normalized_expression_filtered[het_mOSN_UMI < 2] <- 0
write.csv(het_mOSN_normalized_expression_filtered, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv")
# cell number table
cellNumber_OR_het_mOSN <- as.data.frame(table(colSums(het_mOSN_normalized_expression_filtered > 0)))
colnames(cellNumber_OR_het_mOSN) <- c("number", "het_mOSN")

### mOSN_homo
homo_mOSN <- subset(x = trim66_harmony, idents = "mOSN_homo")
homo_mOSN_functional_Olfr_Taar_expression <- homo_mOSN[OR_and_taar_functional_reference, ] # 1036 
# normalized expression matrix
homo_mOSN_normalized_expression <- GetAssayData(object = homo_mOSN_functional_Olfr_Taar_expression, slot = "data") #(eg. “counts”, “data”, or “scale.data”). 
homo_mOSN_normalized_expression <- as.matrix(homo_mOSN_normalized_expression)
# UMI matrix
homo_mOSN_UMI <- GetAssayData(object = homo_mOSN_functional_Olfr_Taar_expression, slot = "counts") #(eg. “counts”, “data”, or “scale.data”). 
homo_mOSN_UMI <- as.matrix(homo_mOSN_UMI)
all(colnames(homo_mOSN_normalized_expression) == colnames(homo_mOSN_UMI))
# filter (UMI > 1)
homo_mOSN_normalized_expression_filtered <-homo_mOSN_normalized_expression
homo_mOSN_normalized_expression_filtered[homo_mOSN_UMI < 2] <- 0
write.csv(homo_mOSN_normalized_expression_filtered, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv")
# cell number table
cellNumber_OR_homo_mOSN <- as.data.frame(table(colSums(homo_mOSN_normalized_expression_filtered > 0)))
colnames(cellNumber_OR_homo_mOSN) <- c("number", "homo_mOSN")

# merge cell number tables
cellNumber_OR_mOSN <- full_join(cellNumber_OR_het_mOSN, cellNumber_OR_homo_mOSN)
cellNumber_OR_mOSN[is.na(cellNumber_OR_mOSN)] <- 0
write.csv(cellNumber_OR_mOSN, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/cellNumber_functional_OR_mOSN_filtered_by_UMI_2.csv")

### het_mOSN_single
row <- as.data.frame(t(colSums(het_mOSN_normalized_expression_filtered>0)))
rownames(row) <- c("count")
het_mOSN_normalized_expression_filtered_single <- het_mOSN_normalized_expression_filtered_single[, (het_mOSN_normalized_expression_filtered_single["count", ] == 1)]
het_mOSN_normalized_expression_filtered_single <- het_mOSN_normalized_expression_filtered_single[-1037, ]
write.csv(het_mOSN_normalized_expression_filtered_single, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv")
# cellNumber, het_mOSN_mean, percent
single_het_mOSN_normalized_expression_filtered_level <- as.data.frame(apply(het_mOSN_normalized_expression_filtered_single, 1, sum))
colnames(single_het_mOSN_normalized_expression_filtered_level) <- c("het_mOSN_sum")
single_het_mOSN_normalized_expression_filtered_level$cellNumber <- rowSums(het_mOSN_normalized_expression_filtered_single>0)
single_het_mOSN_normalized_expression_filtered_level$het_mOSN_mean <- single_het_mOSN_normalized_expression_filtered_level$het_mOSN_sum/single_het_mOSN_normalized_expression_filtered_level$cellNumber
single_het_mOSN_normalized_expression_filtered_level$percent <- single_het_mOSN_normalized_expression_filtered_level$cellNumber/ncol(het_mOSN_normalized_expression_filtered_single)
sum(single_het_mOSN_normalized_expression_filtered_level$percent)
write.csv(single_het_mOSN_normalized_expression_filtered_level, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_single_functional_Taar_OR_expression_level.csv")

### homo_mOSN_single
row <- as.data.frame(t(colSums(homo_mOSN_normalized_expression_filtered>0)))
rownames(row) <- c("count")
homo_mOSN_normalized_expression_filtered_single <- rbind(homo_mOSN_normalized_expression_filtered, row)
homo_mOSN_normalized_expression_filtered_single <- homo_mOSN_normalized_expression_filtered_single[, (homo_mOSN_normalized_expression_filtered_single["count", ] == 1)]
homo_mOSN_normalized_expression_filtered_single <- homo_mOSN_normalized_expression_filtered_single[-1037, ]
write.csv(homo_mOSN_normalized_expression_filtered_single, file = "single_functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv")
# cellNumber, homo_mOSN_mean, percent
single_homo_mOSN_normalized_expression_filtered_level <- as.data.frame(apply(homo_mOSN_normalized_expression_filtered_single, 1, sum))
colnames(single_homo_mOSN_normalized_expression_filtered_level) <- c("homo_mOSN_sum")
single_homo_mOSN_normalized_expression_filtered_level$cellNumber <- rowSums(homo_mOSN_normalized_expression_filtered_single>0)
single_homo_mOSN_normalized_expression_filtered_level$homo_mOSN_mean <- single_homo_mOSN_normalized_expression_filtered_level$homo_mOSN_sum/single_homo_mOSN_normalized_expression_filtered_level$cellNumber
single_homo_mOSN_normalized_expression_filtered_level$percent <- single_homo_mOSN_normalized_expression_filtered_level$cellNumber/ncol(homo_mOSN_normalized_expression_filtered_single)
sum(single_homo_mOSN_normalized_expression_filtered_level$percent)
write.csv(single_homo_mOSN_normalized_expression_filtered_level, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_mOSN_single_functional_Taar_OR_expression_level.csv")

### het_mOSN_multiple
row <- as.data.frame(t(colSums(het_mOSN_normalized_expression_filtered>0)))
rownames(row) <- c("count")
het_mOSN_normalized_expression_filtered_multiple <- rbind(het_mOSN_normalized_expression_filtered, row)
het_mOSN_normalized_expression_filtered_multiple <- het_mOSN_normalized_expression_filtered_multiple[, (het_mOSN_normalized_expression_filtered_multiple["count", ]>1)]
het_mOSN_normalized_expression_filtered_multiple <- het_mOSN_normalized_expression_filtered_multiple[-1037, ]
write.csv(het_mOSN_normalized_expression_filtered_multiple, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/multiple_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv")
# cellNumber, het_mOSN_mean, percent
multiple_het_mOSN_normalized_expression_filtered_level <- as.data.frame(apply(het_mOSN_normalized_expression_filtered_multiple, 1, sum))
colnames(multiple_het_mOSN_normalized_expression_filtered_level) <- c("het_mOSN_sum")
multiple_het_mOSN_normalized_expression_filtered_level$cellNumber <- rowSums(het_mOSN_normalized_expression_filtered_multiple>0)
multiple_het_mOSN_normalized_expression_filtered_level$het_mOSN_mean <- multiple_het_mOSN_normalized_expression_filtered_level$het_mOSN_sum/multiple_het_mOSN_normalized_expression_filtered_level$cellNumber
multiple_het_mOSN_normalized_expression_filtered_level$percent <- multiple_het_mOSN_normalized_expression_filtered_level$cellNumber/ncol(het_mOSN_normalized_expression_filtered_multiple)
sum(multiple_het_mOSN_normalized_expression_filtered_level$percent) # 2.103333
write.csv(multiple_het_mOSN_normalized_expression_filtered_level, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_multiple_functional_Taar_OR_expression_level.csv")

### homo_mOSN_multiple
row <- as.data.frame(t(colSums(homo_mOSN_normalized_expression_filtered>0)))
rownames(row) <- c("count")
homo_mOSN_normalized_expression_filtered_multiple <- rbind(homo_mOSN_normalized_expression_filtered, row)
homo_mOSN_normalized_expression_filtered_multiple <- homo_mOSN_normalized_expression_filtered_multiple[, (homo_mOSN_normalized_expression_filtered_multiple["count", ]>1)]
homo_mOSN_normalized_expression_filtered_multiple <- homo_mOSN_normalized_expression_filtered_multiple[-1037, ]
write.csv(homo_mOSN_normalized_expression_filtered_multiple, file = "multiple_functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv")
# cellNumber, homo_mOSN_mean, percent
multiple_homo_mOSN_normalized_expression_filtered_level <- as.data.frame(apply(homo_mOSN_normalized_expression_filtered_multiple, 1, sum))
colnames(multiple_homo_mOSN_normalized_expression_filtered_level) <- c("homo_mOSN_sum")
multiple_homo_mOSN_normalized_expression_filtered_level$cellNumber <- rowSums(homo_mOSN_normalized_expression_filtered_multiple>0)
multiple_homo_mOSN_normalized_expression_filtered_level$homo_mOSN_mean <- multiple_homo_mOSN_normalized_expression_filtered_level$homo_mOSN_sum/multiple_homo_mOSN_normalized_expression_filtered_level$cellNumber
multiple_homo_mOSN_normalized_expression_filtered_level$percent <- multiple_homo_mOSN_normalized_expression_filtered_level$cellNumber/ncol(homo_mOSN_normalized_expression_filtered_multiple)
sum(multiple_homo_mOSN_normalized_expression_filtered_level$percent) # 7.051103
write.csv(multiple_homo_mOSN_normalized_expression_filtered_level, file = "homo_mOSN_multiple_functional_Taar_OR_expression_level.csv")

##### step 7.2: imOSN (UMI > 0) het/homo single/multiple 受体表达矩阵 (only_functional_receptors) --------------------------------------------------------------------
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"
### het_imOSN
het_imOSN <- subset(x = trim66_harmony, idents = "imOSN_het")
het_imOSN_OR_expression <- het_imOSN[OR_and_taar_functional_reference, ] # 1036 
dim(het_imOSN_OR_expression)
# normalized expression matrix
het_imOSN_normalized_expression <- GetAssayData(object = het_imOSN_OR_expression, slot = "data") #(eg. “counts”, “data”, or “scale.data”). 
het_imOSN_normalized_expression <- as.matrix(het_imOSN_normalized_expression)
het_imOSN_normalized_expression_filtered <-het_imOSN_normalized_expression[colSums(het_imOSN_normalized_expression)> 0, ]
# cell number
cellNumber_het <- as.data.frame(t(colSums(het_imOSN_normalized_expression>0)))
cellNumber_het_df <- as.data.frame(table(as.numeric(cellNumber_het)))
colnames(cellNumber_het_df) <- c("number", "het")

### imOSN_homo
homo_imOSN <- subset(x = trim66_harmony, idents = "imOSN_homo")
homo_imOSN_OR_expression <- homo_imOSN[OR_and_taar_functional_reference, ] # 1036 
# normalized expression matrix
homo_imOSN_normalized_expression <- GetAssayData(object = homo_imOSN_OR_expression, slot = "data") #(eg. “counts”, “data”, or “scale.data”). 
homo_imOSN_normalized_expression <- as.matrix(homo_imOSN_normalized_expression)
cellNumber_homo <- as.data.frame(t(colSums(homo_imOSN_normalized_expression>0)))
rownames(cellNumber_homo) <- c("count")
homo_imOSN_normalized_expression_filtered <- rbind(homo_imOSN_normalized_expression, cellNumber_homo)
homo_imOSN_normalized_expression_filtered <- homo_imOSN_normalized_expression_filtered[, (homo_imOSN_normalized_expression_filtered["count", ]>0)]
homo_imOSN_normalized_expression_filtered <- homo_imOSN_normalized_expression_filtered[-1037, ]
# cell number
cellNumber_homo_df <- as.data.frame(table(as.numeric(cellNumber_homo)))
colnames(cellNumber_homo_df) <- c("number", "homo")

# merge cell number tables
cellNumber_OR_imOSN <- full_join(cellNumber_het_df, cellNumber_homo_df)
cellNumber_OR_imOSN[is.na(cellNumber_OR_imOSN)] <- 0

write.csv(cellNumber_OR_imOSN, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/cellNumber_functional_OR_imOSN.csv")
write.csv(het_imOSN_normalized_expression_filtered, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_imOSN_functional_normalized_expression.csv")
write.csv(homo_imOSN_normalized_expression_filtered, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_imOSN_functional_normalized_expression.csv")

### het_imOSN_single
row <- as.data.frame(t(colSums(het_imOSN_normalized_expression>0)))
rownames(row) <- c("count")
het_imOSN_normalized_expression_single <- rbind(het_imOSN_normalized_expression, row)
het_imOSN_normalized_expression_single <- het_imOSN_normalized_expression_single[, (het_imOSN_normalized_expression_single["count", ] == 1)]
het_imOSN_normalized_expression_single <- het_imOSN_normalized_expression_single[-nrow(het_imOSN_normalized_expression_single), ]
write.csv(het_imOSN_normalized_expression_single, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/single_functional_or_het_imOSN_normalized_expression.csv")

### homo_imOSN_single
row <- as.data.frame(t(colSums(homo_imOSN_normalized_expression>0)))
rownames(row) <- c("count")
homo_imOSN_normalized_expression_single <- rbind(homo_imOSN_normalized_expression, row)
homo_imOSN_normalized_expression_single <- homo_imOSN_normalized_expression_single[, (homo_imOSN_normalized_expression_single["count", ] == 1)]
homo_imOSN_normalized_expression_single <- homo_imOSN_normalized_expression_single[-nrow(homo_imOSN_normalized_expression_single), ]
write.csv(homo_imOSN_normalized_expression_single, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/single_functional_or_homo_imOSN_normalized_expression.csv")

### het_imOSN_multiple
row <- as.data.frame(t(colSums(het_imOSN_normalized_expression>0)))
rownames(row) <- c("count")
het_imOSN_normalized_expression_multiple <- rbind(het_imOSN_normalized_expression, row)
het_imOSN_normalized_expression_multiple <- het_imOSN_normalized_expression_multiple[, (het_imOSN_normalized_expression_multiple["count", ]>1)]
het_imOSN_normalized_expression_multiple <- het_imOSN_normalized_expression_multiple[-nrow(het_imOSN_normalized_expression_multiple), ]
write.csv(het_imOSN_normalized_expression_multiple, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/multiple_functional_or_het_imOSN_normalized_expression.csv")

### homo_imOSN_multiple
row <- as.data.frame(t(colSums(homo_imOSN_normalized_expression>0)))
rownames(row) <- c("count")
homo_imOSN_normalized_expression_multiple <- rbind(homo_imOSN_normalized_expression, row)
homo_imOSN_normalized_expression_multiple <- homo_imOSN_normalized_expression_multiple[, (homo_imOSN_normalized_expression_multiple["count", ]>1)]
homo_imOSN_normalized_expression_multiple <- homo_imOSN_normalized_expression_multiple[-nrow(homo_imOSN_normalized_expression_multiple), ]
write.csv(homo_imOSN_normalized_expression_multiple, file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/multiple_functional_or_homo_imOSN_normalized_expression.csv")

#### Step 7: single/multiple receptor percent (mOSN) ####
##### Plot 7.1: barPlot: the percent of mOSNs expressing a particular receptor (the top 50 in homo_mOSN) (single/multiple) (only_functional_receptors)####
single_het_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_single_functional_Taar_OR_expression_level.csv", header = T, row.names = 1)
single_homo_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_mOSN_single_functional_Taar_OR_expression_level.csv", header = T, row.names = 1)

### single
colnames(single_het_mOSN_normalized_expression_filtered_level) <- c("het_mOSN_sum", "het_mOSN_cellNumber", "het_mOSN_mean", "het_mOSN_percent")
colnames(single_homo_mOSN_normalized_expression_filtered_level) <- c("homo_mOSN_sum", "homo_mOSN_cellNumber", "homo_mOSN_mean", "homo_mOSN_percent")
mOSN_single_Taar_OR_expression_level <- merge(single_het_mOSN_normalized_expression_filtered_level, single_homo_mOSN_normalized_expression_filtered_level, by = 0)
mOSN_single_Taar_OR_expression_level <- mOSN_single_Taar_OR_expression_level[order(mOSN_single_Taar_OR_expression_level$homo_mOSN_percent, decreasing = T), ]
data_plot <- mOSN_single_Taar_OR_expression_level[1:50, ] # homo排列的前50
data_plot[is.na(data_plot)] <- 0

data_plot_het <- data_plot[, c("Row.names", "het_mOSN_percent")]
data_plot_het$class <- c("het")
colnames(data_plot_het) <- c("identity", "percent", "class")

data_plot_homo <- data_plot[, c("Row.names", "homo_mOSN_percent")]
data_plot_homo$class <-c("homo") 
colnames(data_plot_homo) <- c("identity", "percent", 'class') 

data_plot <- rbind(data_plot_homo, data_plot_het)
data_plot_homo <- data_plot_homo[order(data_plot_homo$percent, decreasing = T), ]
data_plot$identity <- factor(data_plot$identity, level = data_plot_homo$identity) #按homo前50排列

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/percent_cell_expressing_functional_TAAR_OR_in_single_mOSN.pdf", width = 20, height = 7)
ggplot(data = data_plot, aes(identity, percent))+
  geom_bar(aes(color = class, fill = class), stat = "identity", position = "dodge")+
  scale_fill_manual(values = c("springgreen3", "lightcoral"))+ #改变柱状图轮廓的颜色
  scale_colour_manual(values = c("springgreen3", "lightcoral"))+ #改变柱状图内填充的颜色
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.y = element_text(size = 18, 
                                    # family = "myFont", 
                                    color = "black", 
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90), 
        axis.title.x = element_text(size = 18, 
                                    # family = "myFont", 
                                    color = "black", 
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 0), 
        legend.title = element_text(color = "black", # 修改图例的标题
                                    size = 16, 
                                    face = "plain"), 
        legend.text = element_text(color = "black", # 设置图例标签文字
                                   size = 14, 
                                   face = "plain"), 
        axis.text.x = element_text(size = 16, # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "plain", # face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), # 角度
        axis.text.y = element_text(size = 16, 
                                   color = "black", 
                                   face = "plain", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0))+
  scale_y_continuous(expand = c(0, 0))
dev.off()

### multiple
multiple_het_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_multiple_Taar_OR_expression_level.csv", header = T, row.names = 1)
multiple_homo_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_mOSN_multiple_Taar_OR_expression_level.csv", header = T, row.names = 1)

colnames(multiple_het_mOSN_normalized_expression_filtered_level) <- c("het_mOSN_sum", "het_mOSN_cellNumber", "het_mOSN_mean", "het_mOSN_percent")
colnames(multiple_homo_mOSN_normalized_expression_filtered_level) <- c("homo_mOSN_sum", "homo_mOSN_cellNumber", "homo_mOSN_mean", "homo_mOSN_percent")
mOSN_multiple_Taar_OR_expression_level <- merge(multiple_het_mOSN_normalized_expression_filtered_level, multiple_homo_mOSN_normalized_expression_filtered_level, by = 0)
mOSN_multiple_Taar_OR_expression_level <- mOSN_multiple_Taar_OR_expression_level[order(mOSN_multiple_Taar_OR_expression_level$homo_mOSN_percent, decreasing = T), ]
data_plot <- mOSN_multiple_Taar_OR_expression_level[1:50, ]
data_plot[is.na(data_plot)] <- 0

data_plot_het <- data_plot[, c("Row.names", "het_mOSN_percent")]
data_plot_het$class <- c("het")
colnames(data_plot_het) <- c("identity", "percent", "class")

data_plot_homo <- data_plot[, c("Row.names", "homo_mOSN_percent")]
data_plot_homo$class <-c("homo") 
colnames(data_plot_homo) <- c("identity", "percent", 'class') 

data_plot <- rbind(data_plot_homo, data_plot_het)
data_plot_homo <- data_plot_homo[order(data_plot_homo$percent, decreasing = T), ]
data_plot$identity <- factor(data_plot$identity, level = data_plot_homo$identity)

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/percent_cell_expressing_functional_TAAR_OR_in_multiple_mOSN.pdf", width = 20, height = 7)
ggplot(data = data_plot, aes(identity, percent))+
  geom_bar(aes(color = class, fill = class), stat = "identity", position = "dodge")+
  scale_fill_manual(values = c("springgreen3", "lightcoral"))+ #改变柱状图轮廓的颜色
  scale_colour_manual(values = c("springgreen3", "lightcoral"))+ #改变柱状图内填充的颜色
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.y = element_text(size = 18, 
                                    # family = "myFont", 
                                    color = "black", 
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90), 
        axis.title.x = element_text(size = 18, 
                                    # family = "myFont", 
                                    color = "black", 
                                    face = "plain", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 0), 
        legend.title = element_text(color = "black", # 修改图例的标题
                                    size = 16, 
                                    face = "plain"), 
        legend.text = element_text(color = "black", # 设置图例标签文字
                                   size = 14, 
                                   face = "plain"), 
        axis.text.x = element_text(size = 16, # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "plain", # face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 90), # 角度
        axis.text.y = element_text(size = 16, 
                                   color = "black", 
                                   face = "plain", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0))+
  scale_y_continuous(expand = c(0, 0))
dev.off()

##### Plot 7.2: circosPlot: the percentage of mOSNs expressing a particular receptor (all receptors) (single) (only_functional_receptors) #####
single_het_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_single_functional_Taar_OR_expression_level.csv", header = T, row.names = 1)
single_homo_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_mOSN_single_functional_Taar_OR_expression_level.csv", header = T, row.names = 1)

### data preparation
#1. extract olfactory receptor information from gencode 
gencode_information <- read.delim("Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, sep = "\t", row.names = NULL, as.is = T)
# extract taar & olfr gene
OR_information_gencode <- gencode_information[grep("^Olfr|^Taar", gencode_information$gene_name), ] # 1428 rows
# only save protein_coding gene
OR_information_gencode_only_functional <- OR_information_gencode[grep("protein_coding", OR_information_gencode$gene_type), ]#1155
# test: OR_information_gencode_only_functional[which(OR_information_gencode_only_functional$gene_name == "Taar7c"), ]

# 2. add genomic_start_relative
OR_information_gencode_only_functional$genomic_start_relative <- 1
for (i in unique(OR_information_gencode_only_functional$chr_name)) {
  OR_information_gencode_only_functional[which(OR_information_gencode_only_functional$chr_name == i), ]$genomic_start_relative <- 1:sum(OR_information_gencode_only_functional$chr_name == i)
}

# 3. add gemomic_end_relative
OR_information_gencode_only_functional$genomic_end_relative <- OR_information_gencode_only_functional$genomic_start_relative + 0.8

# # 4.extract class information from lq All_OR_information
# class_information <- read.csv("All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv", header = T, row.names = 1)
# OR_information_gencode_only_functional$OR_class <- class_information$OR_Class[match(OR_information_gencode_only_functional$gene_name, rownames(class_information))]
# OR_information_gencode_only_functional$OR_class[grep("^Taar", OR_information_gencode_only_functional$gene_name)] <- "Taar"

# 5.export necessary column
OR_information_circos_data <- OR_information_gencode_only_functional[, c("chr_name", "genomic_start_relative", "genomic_end_relative", "gene_name")]

het_homo_mOSN_percent_circos_data <- OR_information_circos_data
het_homo_mOSN_percent_circos_data$het <- single_het_mOSN_normalized_expression_filtered_level[match(het_homo_mOSN_percent_circos_data$gene_name, rownames(single_het_mOSN_normalized_expression_filtered_level)), "percent"]
het_homo_mOSN_percent_circos_data$homo <- single_homo_mOSN_normalized_expression_filtered_level[match(het_homo_mOSN_percent_circos_data$gene_name, rownames(single_homo_mOSN_normalized_expression_filtered_level)), "percent"]
het_homo_mOSN_percent_circos_data[is.na(het_homo_mOSN_percent_circos_data)] <- 0
het_homo_mOSN_percent_circos_data <- het_homo_mOSN_percent_circos_data[-c(1154, 1155), ]

### plot
# color_brewer
# http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
# https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21) 
col_for_chr_percent <- col_for_chr_universal
# obtain the colors for chr used in step1.1
chr_all <- paste0("chr", c(1:19, "X", "Y"))
col_for_chr_percent <- col_for_chr_percent[match(chr_all, unique(het_homo_mOSN_percent_circos_data$chr_name))]
col_for_chr_percent <- col_for_chr_percent[!is.na(col_for_chr_percent)]

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/single_percent_circosPlot_percent_het_homo_mOSN_with_top20.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(het_homo_mOSN_percent_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels
# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_percent, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1

# 2. add receptor name
homo_top20_receptor <- top_n(het_homo_mOSN_percent_circos_data[order(het_homo_mOSN_percent_circos_data$homo, decreasing = T), ], 20)
circos.genomicLabels(homo_top20_receptor, labels.column = 4, side = "outside")

# 3. add histgram
# add the second track of real data
circos.genomicTrack(het_homo_mOSN_percent_circos_data, 
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop.column = 2, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original data.
                                         border = "springgreen3", # otherwise there are black border
                                         col = "springgreen3") #use [[]] instead of [] to obtain vector other than dataframe.
                      circos.genomicRect(region, value, 
                                         ytop.column = 3, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original data.
                                         border = "lightcoral", # otherwise there are black border
                                         col = "lightcoral") #use [[]] instead of [] to obtain vector other than dataframe.
                    })
circos.clear()
dev.off()

##### Plot 7.3: circosPlot: the percentage of mOSNs expressing a particular receptor (all receptors) (multiple) (only_functional_receptors) #####
multiple_het_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_multiple_functional_Taar_OR_expression_level.csv", header = T, row.names = 1)
multiple_homo_mOSN_normalized_expression_filtered_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_mOSN_multiple_functional_Taar_OR_expression_level.csv", header = T, row.names = 1)
### data preparation
#1. extract olfactory receptor information from gencode 
gencode_information <- read.delim("Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, sep = "\t", row.names = NULL, as.is = T)
# extract taar & olfr gene
OR_information_gencode <- gencode_information[grep("^Olfr|^Taar", gencode_information$gene_name), ] # 1428 rows
# only save protein_coding gene
OR_information_gencode_only_functional <- OR_information_gencode[grep("protein_coding", OR_information_gencode$gene_type), ]#1155
# test: OR_information_gencode_only_functional[which(OR_information_gencode_only_functional$gene_name == "Taar7c"), ]

# 2. add genomic_start_relative
OR_information_gencode_only_functional$genomic_start_relative <- 1
for (i in unique(OR_information_gencode_only_functional$chr_name)) {
  OR_information_gencode_only_functional[which(OR_information_gencode_only_functional$chr_name == i), ]$genomic_start_relative <- 1:sum(OR_information_gencode_only_functional$chr_name == i)
}

# 3. add gemomic_end_relative
OR_information_gencode_only_functional$genomic_end_relative <- OR_information_gencode_only_functional$genomic_start_relative + 0.8

# # 4.extract class information from lq All_OR_information
# class_information <- read.csv("All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv", header = T, row.names = 1)
# OR_information_gencode_only_functional$OR_class <- class_information$OR_Class[match(OR_information_gencode_only_functional$gene_name, rownames(class_information))]
# OR_information_gencode_only_functional$OR_class[grep("^Taar", OR_information_gencode_only_functional$gene_name)] <- "Taar"

# 5.export necessary column
OR_information_circos_data <- OR_information_gencode_only_functional[, c("chr_name", "genomic_start_relative", "genomic_end_relative", "gene_name")]

het_homo_mOSN_percent_circos_data <- OR_information_circos_data
het_homo_mOSN_percent_circos_data$het <- multiple_het_mOSN_normalized_expression_filtered_level[match(het_homo_mOSN_percent_circos_data$gene_name, rownames(multiple_het_mOSN_normalized_expression_filtered_level)), "percent"]
het_homo_mOSN_percent_circos_data$homo <- multiple_homo_mOSN_normalized_expression_filtered_level[match(het_homo_mOSN_percent_circos_data$gene_name, rownames(multiple_homo_mOSN_normalized_expression_filtered_level)), "percent"]
het_homo_mOSN_percent_circos_data[is.na(het_homo_mOSN_percent_circos_data)] <- 0
het_homo_mOSN_percent_circos_data <- het_homo_mOSN_percent_circos_data[-c(1154, 1155), ]

### plot
# color_brewer
# http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
# https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21) 
col_for_chr_percent <- col_for_chr_universal
# obtain the colors for chr used in step1.1
chr_all <- paste0("chr", c(1:19, "X", "Y"))
col_for_chr_percent <- col_for_chr_percent[match(chr_all, unique(het_homo_mOSN_percent_circos_data$chr_name))]
col_for_chr_percent <- col_for_chr_percent[!is.na(col_for_chr_percent)]

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/multiple_percent_circosPlot_percent_het_homo_mOSN_with_top20.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(het_homo_mOSN_percent_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels
# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_percent, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1

# 2. add receptor name
homo_top20_receptor <- top_n(het_homo_mOSN_percent_circos_data[order(het_homo_mOSN_percent_circos_data$homo, decreasing = T), ], 20)
circos.genomicLabels(homo_top20_receptor, labels.column = 4, side = "outside")

# 3. add histgram
# add the second track of real data
circos.genomicTrack(het_homo_mOSN_percent_circos_data, 
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop.column = 2, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original data.
                                         border = "springgreen3", # otherwise there are black border
                                         col = "springgreen3") #use [[]] instead of [] to obtain vector other than dataframe.
                      circos.genomicRect(region, value, 
                                         ytop.column = 3, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original data.
                                         border = "lightcoral", # otherwise there are black border
                                         col = "lightcoral") #use [[]] instead of [] to obtain vector other than dataframe.
                    })
circos.clear()
dev.off()

#### step 8: 0/single/multiple het/homo mOSN/imOSN ####
##### plot 8.1: stack bar plot of the percent of mOSNs expressing 0/single/multiple receptors (het/homo) (only_functional_receptors) ####
cellNumber_OR_mOSN <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/cellNumber_functional_OR_mOSN_filtered_by_UMI_2.csv", header = T, row.names = 1)
rownames(cellNumber_OR_mOSN) <- cellNumber_OR_mOSN$number
cellNumber_OR_mOSN$het_mOSN_percent <- cellNumber_OR_mOSN$het/sum(cellNumber_OR_mOSN$het)
cellNumber_OR_mOSN$homo_mOSN_percent <- cellNumber_OR_mOSN$homo/sum(cellNumber_OR_mOSN$homo)

cellnumber_new <- cellNumber_OR_mOSN[1:2, ]
cellnumber_new["multiple", ] <- colSums(cellNumber_OR_mOSN[3:nrow(cellNumber_OR_mOSN), ])
fisher.test(cellnumber_new[2:3, 2:3])
fisher.test(cellnumber_new[2:3, 2:3])$p.value

barplot_data <- cellnumber_new[1:3, 2:3] %>% 
  rownames_to_column(var = 'number') %>% 
  pivot_longer( cols = c("het", "homo"), 
                names_to = 'sampletype', 
                values_to = 'cellCount')


pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/stack_barPlot_mOSN_functional_receptor_number.pdf", width = 5, height = 6)

p1 <- ggplot(barplot_data[which(barplot_data$sampletype == "het"), ], aes(x = sampletype, y = cellCount, fill = number)) +
  geom_bar(position = 'stack', stat = 'identity', color = "black", fill = c("azure4", "deepskyblue2", "hotpink2"))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  xlab(NULL)+
  theme (axis.line.x = element_blank (), 
         axis.ticks.x = element_blank (), 
         axis.line.y = element_blank (), 
         axis.text.y = element_blank (), # remove y axis labels 
         axis.ticks.y = element_blank () )+ #remove y axis ticks 
  ylab(NULL)


p2 <- p1 + ggplot(barplot_data[which(barplot_data$sampletype == "homo"), ], aes(x = sampletype, y = cellCount, fill = number)) +
  geom_bar(position = 'stack', stat = 'identity', color = "black", fill = c("azure4", "deepskyblue2", "hotpink2"))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  xlab(NULL)+
  theme (axis.line.x = element_blank (), 
         axis.ticks.x = element_blank (), 
         axis.line.y = element_blank (), 
         axis.text.y = element_blank (), # remove y axis labels 
         axis.ticks.y = element_blank () )+ #remove y axis ticks 
  ylab(NULL)

print(p2)
dev.off()
##### plot 8.2: stack bar plot of the percent of imOSNs expressing 0/single/multiple receptors (het/homo) (only_functional_receptors) ####
cellNumber_OR_imOSN <- read.csv("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/cellNumber_functional_OR_imOSN.csv", header = T, row.names = 1)
rownames(cellNumber_OR_imOSN) <- cellNumber_OR_imOSN$number
cellNumber_OR_imOSN$het_imOSN_percent <- cellNumber_OR_imOSN$het/sum(cellNumber_OR_imOSN$het)
cellNumber_OR_imOSN$homo_imOSN_percent <- cellNumber_OR_imOSN$homo/sum(cellNumber_OR_imOSN$homo)

cellnumber_new <- cellNumber_OR_imOSN[1:2, ]
cellnumber_new["multiple", ] <- colSums(cellNumber_OR_imOSN[3:nrow(cellNumber_OR_imOSN), ])
fisher.test(cellnumber_new[2:3, 2:3])
fisher.test(cellnumber_new[2:3, 2:3])$p.value

# Fisher's Exact Test for Count Data

# data: cellnumber_new[2:3, 2:3]
# p-value = 2.505e-05
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
# 1.326896 2.178459
# sample estimates:
# odds ratio 
# 1.704976 


barplot_data <- cellnumber_new[1:3, 2:3] %>% 
  rownames_to_column(var = 'number') %>% 
  pivot_longer( cols = c("het", "homo"), 
                names_to = 'sampletype', 
                values_to = 'cellCount')

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/stack_barPlot_imOSN_functional_receptor_number.pdf", width = 5, height = 6)

p1 <- ggplot(barplot_data[which(barplot_data$sampletype == "het"), ], aes(x = sampletype, y = cellCount, fill = number)) +
  geom_bar(position = 'stack', stat = 'identity', color = "black", fill = c("azure4", "deepskyblue2", "hotpink2"))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  xlab(NULL)+
  theme (axis.line.x = element_blank (), 
         axis.ticks.x = element_blank (), 
         axis.line.y = element_blank (), 
         axis.text.y = element_blank (), # remove y axis labels 
         axis.ticks.y = element_blank () )+ #remove y axis ticks 
  ylab(NULL)


p2 <- p1 + ggplot(barplot_data[which(barplot_data$sampletype == "homo"), ], aes(x = sampletype, y = cellCount, fill = number)) +
  geom_bar(position = 'stack', stat = 'identity', color = "black", fill = c("azure4", "deepskyblue2", "hotpink2"))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  xlab(NULL)+
  theme (axis.line.x = element_blank (), 
         axis.ticks.x = element_blank (), 
         axis.line.y = element_blank (), 
         axis.text.y = element_blank (), # remove y axis labels 
         axis.ticks.y = element_blank () )+ #remove y axis ticks 
  ylab(NULL)

print(p2)
dev.off()

#### step 9: mOSN/imOSN het/homo UMI violin plot##########################################
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"

het_mOSN <- subset(x = trim66_harmony, idents = "mOSN_het")
het_mOSN_UMI <- GetAssayData(object = het_mOSN, slot = "counts") # (eg. “counts”, “data”, or “scale.data”)
het_mOSN_UMI <- data.frame(UMI_sum = colSums(het_mOSN_UMI)) 
het_mOSN_UMI$genotype_celltype <- "het_mOSN_UMI_sum"

homo_mOSN <- subset(x = trim66_harmony, idents = "mOSN_homo")
homo_mOSN_UMI <- GetAssayData(object = homo_mOSN, slot = "counts")
homo_mOSN_UMI <- data.frame(UMI_sum = colSums(homo_mOSN_UMI)) 
homo_mOSN_UMI$genotype_celltype <- "homo_mOSN_UMI_sum"

het_imOSN <- subset(x = trim66_harmony, idents = "imOSN_het")
het_imOSN_UMI <- GetAssayData(object = het_imOSN, slot = "counts") 
het_imOSN_UMI <- data.frame(UMI_sum = colSums(het_imOSN_UMI)) 
het_imOSN_UMI$genotype_celltype <- "het_imOSN_UMI_sum"

homo_imOSN <- subset(x = trim66_harmony, idents = "imOSN_homo")
homo_imOSN_UMI <- GetAssayData(object = homo_imOSN, slot = "counts") 
homo_imOSN_UMI <- data.frame(UMI_sum = colSums(homo_imOSN_UMI)) 
homo_imOSN_UMI$genotype_celltype <- "homo_imOSN_UMI_sum"

multirbind <-function(dat = list(), ...){
  if(length(dat) < 2)return(as.data.frame(dat))
  rbinddat <-dat[[1]]
  dat[[1]] <-NULL
  for(i in dat){
    rbinddat <-rbind(rbinddat, i)
  }
  return(rbinddat)
}
UMI_sum_list <- list(het_mOSN_UMI, homo_mOSN_UMI, het_imOSN_UMI, homo_imOSN_UMI)
UMI_sum <- multirbind(UMI_sum_list)
UMI_sum$genotype_celltype <- factor(UMI_sum$genotype_celltype, level = c("het_mOSN_UMI_sum", "homo_mOSN_UMI_sum", "het_imOSN_UMI_sum", "homo_imOSN_UMI_sum"))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/UMI_sum_violinPlot.pdf", width = 10, height = 8)
ggplot(UMI_sum, aes(x = genotype_celltype, y = UMI_sum, fill = genotype_celltype))+
  geom_jitter(size = 0.5)+
  geom_violin()+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  theme_bw()+theme(panel.border = element_blank(), 
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                   axis.line = element_line(colour = "black"))
dev.off()

#### step 10: the distribution of cells which expressing similar receptor #####################
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
# OR info
all_gene_information <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", sep = "\t", header = T)
OR_and_taar_functional_reference <- all_gene_information$gene_name[intersect(grep(pattern = "protein_coding", all_gene_information$gene_type), grep(pattern = "^Olfr|^Taar", all_gene_information$gene_name))]
OR_and_taar_functional_reference <- OR_and_taar_functional_reference[!duplicated(OR_and_taar_functional_reference)] # 1152
OR_and_taar_functional_reference <- intersect(OR_and_taar_functional_reference, rownames(trim66_harmony))

# all
list_olfr_umap <- list()
pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/featurePlot_functional_olfr_umap.pdf", width = 7, height = 7)
for (i in c(1:1036))
{list_olfr_umap[[i]] <- FeaturePlot(object = trim66_harmony, 
                                    features = OR_and_taar_functional_reference[i], 
                                    reduction = "umap")
}
print(list_olfr_umap)
dev.off()

# split by genotype
Idents(trim66_harmony) <- "genotype"
list_olfr_umap <- list()
pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/featurePlot_functional_olfr_umap_het_homo.pdf", width = 14, height = 7)
for (i in c(1:1036)){
  list_olfr_umap[[i]] <- FeaturePlot(object = trim66_harmony, 
                                    features = OR_and_taar_functional_reference[i], 
                                    split.by = "genotype",
                                    reduction = "umap",
                                    order = T)
}
print(list_olfr_umap)
dev.off()

#### step 11: pseudotime analysis ####
##### step 11.1: pseudotime analysis (monocle3) ####
# c("mOSN", "imOSN", "INP", "GBC")
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype"
trim66_harmony_filter <- subset(trim66_harmony, idents = c("mOSN", "imOSN", "INP", "GBC"))
Idents(trim66_harmony_filter) <- "genotype"
het <- subset(trim66_harmony_filter, idents = "het")
homo <- subset(trim66_harmony_filter, idents = "homo")

#### het 
### 创建cds对象 
#提取counts矩阵
het_counts <- GetAssayData(het, assay = 'RNA', slot = 'counts')
#提取metadata
het_metadata <- het@meta.data
#提取gene
gene_annotation <- data.frame(gene_short_name = rownames(het))
rownames(gene_annotation) <- rownames(het)
#创建cds
het_cds <- new_cell_data_set(het_counts, 
                             cell_metadata = het_metadata, 
                             gene_metadata = gene_annotation)
### 归一化处理 
het_cds <- preprocess_cds(het_cds, num_dim = 100)
plot_pc_variance_explained(het_cds)

### 降维 
het_cds <- reduce_dimension(het_cds, preprocess_method = "PCA") #默认使用UMAP 

plot_cells(het_cds, reduction_method = "UMAP", color_cells_by = "celltype", show_trajectory_graph = FALSE) + 
  ggtitle('cds.umap')

#从seurat导入整合过的umap坐标
cds.embed <- het_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(het, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
het_cds@int_colData$reducedDims$UMAP <- int.embed

plot_cells(het_cds, reduction_method = "UMAP", color_cells_by = "celltype", show_trajectory_graph = FALSE) + 
  ggtitle('int.umap')

### cluster 
het_cds <- cluster_cells(het_cds)

p1 <- plot_cells(het_cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(het_cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)

### 拟时序分析 
het_cds <- learn_graph(het_cds, )

plot_cells(het_cds, 
           color_cells_by = "celltype", 
           label_groups_by_cluster = FALSE, 
           label_leaves = FALSE, 
           label_branch_points = FALSE, 
           group_label_size = 4, cell_size = 1.5)
#选择发育起点
het_cds = order_cells(het_cds)

#自动选择发育起点
get_earliest_principal_node <- function(cds, time_bin = "GBC"){
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids, ]))))]
  root_pr_nodes
}

het_cds = order_cells(het_cds, root_pr_nodes = get_earliest_principal_node(het_cds))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/het_pseudotime_1.pdf", width = 16, height = 7)
p1 <- plot_cells(het_cds, 
                 color_cells_by = "celltype", 
                 show_trajectory_graph = FALSE, 
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE, 
                 group_label_size = 4, cell_size = 0.8)+
  xlim(-8, 10)+
  ylim(-11, 8)+
  scale_color_manual(values = c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2"))
p2 <- plot_cells(het_cds, color_cells_by = "pseudotime", 
                 show_trajectory_graph = FALSE, 
                 label_cell_groups = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE, 
                 cell_size = 0.8)+
  xlim(-8, 10)+
  ylim(-11, 8)
grid.arrange(p1, p2, nrow = 1, ncol = 2)
dev.off()

#### homo 
### 创建cds对象 
#提取counts矩阵
homo_counts <- GetAssayData(homo, assay = 'RNA', slot = 'counts')
#提取metadata
homo_metadata <- homo@meta.data
#提取gene
gene_annotation <- data.frame(gene_short_name = rownames(homo))
rownames(gene_annotation) <- rownames(homo)
#创建cds
homo_cds <- new_cell_data_set(homo_counts, 
                              cell_metadata = homo_metadata, 
                              gene_metadata = gene_annotation)
### 归一化处理 
homo_cds <- preprocess_cds(homo_cds, num_dim = 100)
plot_pc_variance_explained(homo_cds)

### 降维 
homo_cds <- reduce_dimension(homo_cds, preprocess_method = "PCA") #默认使用UMAP 
plot_cells(homo_cds, reduction_method = "UMAP", color_cells_by = "celltype", show_trajectory_graph = FALSE) + 
  ggtitle('cds.umap')

#从seurat导入整合过的umap坐标
cds.embed <- homo_cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(homo, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed), ]
homo_cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(homo_cds, reduction_method = "UMAP", color_cells_by = "celltype", show_trajectory_graph = FALSE) + 
  ggtitle('int.umap')

### cluster 
homo_cds <- cluster_cells(homo_cds)
p1 <- plot_cells(homo_cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(homo_cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)

### 拟时序分析 
homo_cds <- learn_graph(homo_cds, )
plot_cells(homo_cds, 
           color_cells_by = "celltype", 
           label_groups_by_cluster = FALSE, 
           label_leaves = FALSE, 
           label_branch_points = FALSE, 
           group_label_size = 4, cell_size = 1.5)
#选择发育起点
homo_cds = order_cells(homo_cds)

#自动选择发育起点
homo_cds = order_cells(homo_cds, root_pr_nodes = get_earliest_principal_node(homo_cds))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/homo_pseudotime_1.pdf", width = 16, height = 7)
p1 <- plot_cells(homo_cds, 
                 color_cells_by = "celltype", 
                 show_trajectory_graph = FALSE, 
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE, 
                 group_label_size = 4, 
                 cell_size = 0.8)+
  xlim(-8, 10)+
  ylim(-11, 8)+
  scale_color_manual(values = c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2"))
p2 <- plot_cells(homo_cds, color_cells_by = "pseudotime", 
                 show_trajectory_graph = FALSE, 
                 label_cell_groups = FALSE, 
                 label_leaves = FALSE, 
                 label_branch_points = FALSE, 
                 cell_size = 0.8)+
  xlim(-8, 10)+
  ylim(-11, 8)
grid.arrange(p1, p2, nrow = 1, ncol = 2)
dev.off()

##### step 11.2.1: pseudotime analysis (slingshot) -----------------------------------
### choose particular celltypes to run 
Idents(trim66_harmony) <- "celltype"
trim66_harmony_filter <- subset(trim66_harmony, idents = c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a"))
Idents(trim66_harmony_filter) <- "genotype"
het_harmony <- subset(trim66_harmony_filter, idents = "het")
homo_harmony <- subset(trim66_harmony_filter, idents = "homo")

### convert objects to SingleCellExperiment objects 
het_sce <- as.SingleCellExperiment(het_harmony, assay = NULL)

### Using Slingshot
het_sce <- slingshot(het_sce, clusterLabels = "celltype", reducedDim = 'UMAP')

# choose reducedDim
reducedDims(het_sce)
# List of length 4
# names(4): PCA HARMONY TSNE UMAP

# choose clustersLabels
colData(het_sce)

summary(het_sce$slingPseudotime_1)

### visuzalize the inferred lineage
col_celltype <- c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2", "cyan3", "palegreen3", "peachpuff3", "lightslateblue")

# without start and end
plot(reducedDims(het_sce)$UMAP, col = col_celltype[het_sce$celltype], pch = 16, asp = 1)
lines(SlingshotDataSet(het_sce), lwd = 2, type = 'lineages', col = 'black')

# with start (HBC) and end (mOSN MV Ms4a)
pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/het_pseudotime_umap.pdf", width = 8, height = 7)
lin2 <- getLineages(reducedDims(het_sce)$UMAP, het_sce$celltype, start.clus = 'HBC', 'GBC', end.clus = c('mOSN', 'Ms4a', 'MV'))
plot(reducedDims(het_sce)$UMAP, col = col_celltype[het_sce$celltype], pch = 16, asp = 1, axes = F, main = "Het")
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
legend("topright", 
       legend = c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a"), 
       fill = col_celltype, # Color of the squares
       border = col_celltype ) # Color of the border of the squares
legend("right", 
       legend = c("start", "end"), 
       fill = c('green', 'red'), # Color of the squares
       border = "black" ) # Color of the border of the squares
dev.off()

### convert objects to SingleCellExperiment objects 
homo_sce <- as.SingleCellExperiment(homo_harmony, assay = NULL)

### Using Slingshot
homo_sce <- slingshot(homo_sce, clusterLabels = "celltype", reducedDim = 'UMAP')

# choose reducedDim
reducedDims(homo_sce)
# List of length 4
# names(4): PCA HARMONY TSNE UMAP

# choose clustersLabels
colData(homo_sce)

summary(homo_sce$slingPseudotime_1)

### visuzalize the inferred lineage

# without start and end 
plot(reducedDims(homo_sce)$UMAP, col = col_celltype, pch = 16, asp = 1)
lines(SlingshotDataSet(homo_sce), lwd = 2, type = 'lineages', col = 'black')

# with start (HBC) and end (mOSN MV Ms4a)
pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/homo_pseudotime_umap.pdf", width = 8, height = 7)
lin2 <- getLineages(reducedDims(homo_sce)$UMAP, homo_sce$celltype, start.clus = 'HBC', 'GBC', end.clus = c('mOSN', 'Ms4a', 'MV'))
plot(reducedDims(homo_sce)$UMAP, col = col_celltype[homo_sce$celltype], pch = 16, asp = 1, axes = F, main = "Homo")
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
legend("topright", 
       legend = c("mOSN", "imOSN", "INP", "GBC", "HBC", "SUS", "MV", "Ms4a"), 
       fill = col_celltype, # Color of the squares
       border = col_celltype ) # Color of the border of the squares
legend("right", 
       legend = c("start", "end"), 
       fill = c('green', 'red'), # Color of the squares
       border = "black" ) # Color of the border of the squares
dev.off()

##### step 11.2.2: pseudotime analysis (slingshot) only with osn lineage ##############
Idents(trim66_harmony) <- "celltype"
trim66_harmony_filter <- subset(trim66_harmony, idents = c("mOSN", "imOSN", "INP", "GBC", "HBC"))
Idents(trim66_harmony_filter) <- "genotype"
het_harmony <- subset(trim66_harmony_filter, idents = "het")
homo_harmony <- subset(trim66_harmony_filter, idents = "homo")

### convert objects to SingleCellExperiment objects 
het_sce <- as.SingleCellExperiment(het_harmony, assay = NULL)

### Using Slingshot
het_sce <- slingshot(het_sce, clusterLabels = "celltype", reducedDim = 'UMAP')

# choose reducedDim
reducedDims(het_sce)
# List of length 4
# names(4): PCA HARMONY TSNE UMAP

# choose clustersLabels
colData(het_sce)

summary(het_sce$slingPseudotime_1)

### visuzalize the inferred lineage
col_celltype <- c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2", "cyan3")

# without start and end
plot(reducedDims(het_sce)$UMAP, col = col_celltype[het_sce$celltype], pch = 16, asp = 1)
lines(SlingshotDataSet(het_sce), lwd = 2, type = 'lineages', col = 'black')

# with start (HBC) and end (mOSN MV Ms4a)
pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/het_pseudotime_umap.pdf", width = 8, height = 7)
lin2 <- getLineages(reducedDims(het_sce)$UMAP, het_sce$celltype, start.clus = 'HBC', end.clus = c('mOSN'))
plot(reducedDims(het_sce)$UMAP, col = col_celltype[het_sce$celltype], pch = 16, asp = 1, axes = F, main = "Het", xlim = c(-8, 10), ylim = c(-11, 8))
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
legend("topright", 
       legend = c("mOSN", "imOSN", "INP", "GBC", "HBC"), 
       fill = col_celltype, # Color of the squares
       border = col_celltype ) # Color of the border of the squares
legend("right", 
       legend = c("start", "end"), 
       fill = c('green', 'red'), # Color of the squares
       border = "black" ) # Color of the border of the squares
dev.off()

### convert objects to SingleCellExperiment objects 
homo_sce <- as.SingleCellExperiment(homo_harmony, assay = NULL)

### Using Slingshot
homo_sce <- slingshot(homo_sce, clusterLabels = "celltype", reducedDim = 'UMAP')

# choose reducedDim
reducedDims(homo_sce)
# List of length 4
# names(4): PCA HARMONY TSNE UMAP

# choose clustersLabels
colData(homo_sce)

summary(homo_sce$slingPseudotime_1)

### visuzalize the inferred lineage
col_celltype <- c("royalblue3", "cornflowerblue", "lightblue3", "paleturquoise2", "cyan3")

# without start and end
plot(reducedDims(homo_sce)$UMAP, col = col_celltype[homo_sce$celltype], pch = 16, asp = 1)
lines(SlingshotDataSet(homo_sce), lwd = 2, type = 'lineages', col = 'black')

# with start (HBC) and end (mOSN MV Ms4a)
pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/homo_pseudotime_umap.pdf", width = 8, height = 7)
lin2 <- getLineages(reducedDims(homo_sce)$UMAP, homo_sce$celltype, start.clus = 'HBC', end.clus = c('mOSN'))
plot(reducedDims(homo_sce)$UMAP, col = col_celltype[homo_sce$celltype], pch = 16, asp = 1, axes = F, main = "homo", xlim = c(-8, 10), ylim = c(-11, 8))+
  xlim(-8, 10)+
  ylim(-11, 8)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)
legend("topright", 
       legend = c("mOSN", "imOSN", "INP", "GBC", "HBC"), 
       fill = col_celltype, # Color of the squares
       border = col_celltype ) # Color of the border of the squares
legend("right", 
       legend = c("start", "end"), 
       fill = c('green', 'red'), # Color of the squares
       border = "black" ) # Color of the border of the squares
dev.off()

#### Step 12: expression violin plot ####
##### plot 12.1: single receptor expression violin plot (single) (only_functional_receptors) ####
### single mOSN
het_mOSN_normalized_expression_filtered_single <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered_single <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/single_functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)

# single expression level
single_expression_level_vlnPlot_data_het_mOSN <- colSums(het_mOSN_normalized_expression_filtered_single)
single_expression_level_vlnPlot_data_homo_mOSN <- colSums(homo_mOSN_normalized_expression_filtered_single)
single_mOSN_list <- list(single_expression_level_vlnPlot_data_het_mOSN, single_expression_level_vlnPlot_data_homo_mOSN)
# Define a function to pad the vectors with NA
pad_vector <- function(x, max_length) {
  c(x, rep(NA, max_length - length(x)))
}
single_mOSN_expression_level_vlnPlot_data <- do.call(cbind, lapply(single_mOSN_list, pad_vector, max_length = max(lengths(single_mOSN_list))))
colnames(single_mOSN_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN")

t.test(single_mOSN_expression_level_vlnPlot_data[, "het_mOSN"], single_mOSN_expression_level_vlnPlot_data[, "homo_mOSN"], paired = FALSE, var.equal = TRUE) 
# t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)

# Two Sample t-test
# 
# data: single_mOSN_expression_level_vlnPlot_data[, "het"] and single_mOSN_expression_level_vlnPlot_data[, "homo"]
# t = 44.516, df = 5700, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# 1.410782 1.540761
# sample estimates:
# mean of x mean of y 
# 4.347045 2.871273 

wilcox.test(single_mOSN_expression_level_vlnPlot_data[, "het_mOSN"], single_mOSN_expression_level_vlnPlot_data[, "homo_mOSN"], paired = F) # 两独立样本的非参数检验

# Wilcoxon rank sum test with continuity correction
# 
# data: single_mOSN_expression_level_vlnPlot_data[, "het"] and single_mOSN_expression_level_vlnPlot_data[, "homo"]
# W = 2955477, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

single_expression_level_vlnPlot_data_longer <- as.data.frame(single_mOSN_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = c("het_mOSN", "homo_mOSN"), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

#### single imOSN
het_imOSN_normalized_expression_single <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/single_functional_or_het_imOSN_normalized_expression.csv", row.names = 1, header = T)
homo_imOSN_normalized_expression_single <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/single_functional_or_homo_imOSN_normalized_expression.csv", row.names = 1, header = T)

#single expression level
single_expression_level_vlnPlot_data_het_imOSN <- colSums(het_imOSN_normalized_expression_single)
single_expression_level_vlnPlot_data_homo_imOSN <- colSums(homo_imOSN_normalized_expression_single)
single_imOSN_list <- list(single_expression_level_vlnPlot_data_het_imOSN, single_expression_level_vlnPlot_data_homo_imOSN)
single_imOSN_expression_level_vlnPlot_data <- do.call(cbind, lapply(single_imOSN_list, pad_vector, max_length = max(lengths(single_imOSN_list))))
colnames(single_imOSN_expression_level_vlnPlot_data) <- c("het_imOSN", "homo_imOSN")

t.test(single_imOSN_expression_level_vlnPlot_data[, "het_imOSN"], single_imOSN_expression_level_vlnPlot_data[, "homo_imOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: single_imOSN_expression_level_vlnPlot_data[, "het"] and single_imOSN_expression_level_vlnPlot_data[, "homo"]
# t = 8.2897, df = 508, p-value = 1.021e-15
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# 0.9651635 1.5647490
# sample estimates:
# mean of x mean of y 
# 3.079684 1.814728

wilcox.test(single_imOSN_expression_level_vlnPlot_data[, "het_imOSN"], single_imOSN_expression_level_vlnPlot_data[, "homo_imOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: single_imOSN_expression_level_vlnPlot_data[, "het"] and single_imOSN_expression_level_vlnPlot_data[, "homo"]
# W = 29532, p-value = 3.417e-10
# alternative hypothesis: true location shift is not equal to 0

single_imOSN_expression_level_vlnPlot_data_longer <- as.data.frame(single_imOSN_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = c("het_imOSN", "homo_imOSN"), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

# single expression level
single_list <- list(single_expression_level_vlnPlot_data_het_mOSN,single_expression_level_vlnPlot_data_homo_mOSN,
               single_expression_level_vlnPlot_data_het_imOSN, single_expression_level_vlnPlot_data_homo_imOSN)

single_expression_level_vlnPlot_data <- do.call(cbind, lapply(single_list, pad_vector, max_length = max(lengths(single_list))))
colnames(single_expression_level_vlnPlot_data) <- c("het_mOSN","homo_mOSN", "het_imOSN", "homo_imOSN")

single_expression_level_vlnPlot_data_longer <- as.data.frame(single_expression_level_vlnPlot_data) %>%  
  pivot_longer( cols =  c("het_mOSN","homo_mOSN", "het_imOSN", "homo_imOSN"),
                names_to = 'sampletype',
                values_to = 'expressionLevel') %>% 
  mutate(sampletype = factor(sampletype, levels = c("het_mOSN","homo_mOSN", "het_imOSN", "homo_imOSN")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/VlnPlot_single_functional_expression_level.pdf",width = 7, height =3)
ggplot(single_expression_level_vlnPlot_data_longer, aes(x=sampletype, y=expressionLevel, fill=sampletype, na.rm = T)) + 
  geom_jitter(size = 0.5)+
  geom_violin(trim = FALSE, na.rm = T)+
  scale_fill_manual(values=alpha(c("springgreen3","lightcoral","springgreen1","lightpink"), 0.8))+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_y_continuous(expand = c(0,0))
dev.off()

##### plot 12.2: top5 receptors expression violin plot (multiple) (only_functional_receptors) ####
het_mOSN_normalized_expression_filtered_multiple <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/multiple_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered_multiple <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/multiple_functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
het_imOSN_normalized_expression <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/multiple_functional_or_het_imOSN_normalized_expression.csv", row.names = 1, header = T)
homo_imOSN_normalized_expression <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/multiple_functional_or_homo_imOSN_normalized_expression.csv", row.names = 1, header = T)

#top1
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
multiple_top1_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(1)(x)])
multiple_top1_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(1)(x)])
top1_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(1)(x)])
top1_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(1)(x)])
top1_list <- list(multiple_top1_expression_level_vlnPlot_data_het_mOSN, multiple_top1_expression_level_vlnPlot_data_homo_mOSN, top1_expression_level_vlnPlot_data_het_imOSN, top1_expression_level_vlnPlot_data_homo_imOSN)
top1_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top1_list, unlist), "length <-", max(lengths(top1_list))))
colnames(top1_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

t.test(top1_expression_level_vlnPlot_data[, "het_mOSN"], top1_expression_level_vlnPlot_data[, "homo_mOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)

# Two Sample t-test
# 
# data: top1_expression_level_vlnPlot_data[, "het_mOSN"] and top1_expression_level_vlnPlot_data[, "homo_mOSN"]
# t = 43.811, df = 9195, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# 1.095917 1.198579
# sample estimates:
# mean of x mean of y 
# 4.292008 3.144760 

wilcox.test(top1_expression_level_vlnPlot_data[, "het_mOSN"], top1_expression_level_vlnPlot_data[, "homo_mOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top1_expression_level_vlnPlot_data[, "het_mOSN"] and top1_expression_level_vlnPlot_data[, "homo_mOSN"]
# W = 6497493, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

t.test(top1_expression_level_vlnPlot_data[, "het_imOSN"], top1_expression_level_vlnPlot_data[, "homo_imOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top1_expression_level_vlnPlot_data[, "het_imOSN"] and top1_expression_level_vlnPlot_data[, "homo_imOSN"]
# t = 11.378, df = 2966, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# 0.6527753 0.9246172
# sample estimates:
# mean of x mean of y 
# 2.890337 2.101641 

wilcox.test(top1_expression_level_vlnPlot_data[, "het_imOSN"], top1_expression_level_vlnPlot_data[, "homo_imOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top1_expression_level_vlnPlot_data[, "het_imOSN"] and top1_expression_level_vlnPlot_data[, "homo_imOSN"]
# W = 646512, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

top1_expression_level_vlnPlot_data_longer <- as.data.frame(top1_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top1_expression_level_vlnPlot_data_longer$sampletype <- factor(top1_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top1_expression_level_vlnPlot_data_longer[top1_expression_level_vlnPlot_data_longer == 0] <- NA

#top2
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
multiple_top2_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(2)(x)])
multiple_top2_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(2)(x)])
top2_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(2)(x)])
top2_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(2)(x)])
top2_list <- list(multiple_top2_expression_level_vlnPlot_data_het_mOSN, multiple_top2_expression_level_vlnPlot_data_homo_mOSN, top2_expression_level_vlnPlot_data_het_imOSN, top2_expression_level_vlnPlot_data_homo_imOSN)
top2_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top2_list, unlist), "length <-", max(lengths(top2_list))))
colnames(top2_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

t.test(top2_expression_level_vlnPlot_data[, "het_mOSN"], top2_expression_level_vlnPlot_data[, "homo_mOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top2_expression_level_vlnPlot_data[, "het_mOSN"] and top2_expression_level_vlnPlot_data[, "homo_mOSN"]
# t = -19.865, df = 9195, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -0.5214020 -0.4277417
# sample estimates:
# mean of x mean of y 
# 1.869583 2.344155
wilcox.test(top2_expression_level_vlnPlot_data[, "het_mOSN"], top2_expression_level_vlnPlot_data[, "homo_mOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top2_expression_level_vlnPlot_data[, "het_mOSN"] and top2_expression_level_vlnPlot_data[, "homo_mOSN"]
# W = 2320819, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
t.test(top2_expression_level_vlnPlot_data[, "het_imOSN"], top2_expression_level_vlnPlot_data[, "homo_imOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top2_expression_level_vlnPlot_data[, "het_imOSN"] and top2_expression_level_vlnPlot_data[, "homo_imOSN"]
# t = -4.5665, df = 2966, p-value = 5.161e-06
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -0.21046032 -0.08401744
# sample estimates:
# mean of x mean of y 
# 0.9073523 1.0545912
wilcox.test(top2_expression_level_vlnPlot_data[, "het_imOSN"], top2_expression_level_vlnPlot_data[, "homo_imOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top2_expression_level_vlnPlot_data[, "het_imOSN"] and top2_expression_level_vlnPlot_data[, "homo_imOSN"]
# W = 372570, p-value = 1.909e-15
# alternative hypothesis: true location shift is not equal to 0

top2_expression_level_vlnPlot_data_longer <- as.data.frame(top2_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top2_expression_level_vlnPlot_data_longer$sampletype <- factor(top2_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top2_expression_level_vlnPlot_data_longer[top2_expression_level_vlnPlot_data_longer == 0] <- NA

#top3
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
multiple_top3_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(3)(x)])
multiple_top3_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(3)(x)])
top3_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(3)(x)])
top3_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(3)(x)])
top3_list <- list(multiple_top3_expression_level_vlnPlot_data_het_mOSN, multiple_top3_expression_level_vlnPlot_data_homo_mOSN, top3_expression_level_vlnPlot_data_het_imOSN, top3_expression_level_vlnPlot_data_homo_imOSN)
top3_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top3_list, unlist), "length <-", max(lengths(top3_list))))
colnames(top3_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

t.test(top3_expression_level_vlnPlot_data[, "het_mOSN"], top3_expression_level_vlnPlot_data[, "homo_mOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top3_expression_level_vlnPlot_data[, "het_mOSN"] and top3_expression_level_vlnPlot_data[, "homo_mOSN"]
# t = -62.072, df = 9195, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -1.719936 -1.614631
# sample estimates:
# mean of x mean of y 
# 0.1276622 1.7949455
wilcox.test(top3_expression_level_vlnPlot_data[, "het_mOSN"], top3_expression_level_vlnPlot_data[, "homo_mOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top3_expression_level_vlnPlot_data[, "het_mOSN"] and top3_expression_level_vlnPlot_data[, "homo_mOSN"]
# W = 541223, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
t.test(top3_expression_level_vlnPlot_data[, "het_imOSN"], top3_expression_level_vlnPlot_data[, "homo_imOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top3_expression_level_vlnPlot_data[, "het_imOSN"] and top3_expression_level_vlnPlot_data[, "homo_imOSN"]
# t = -5.5004, df = 2966, p-value = 4.111e-08
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -0.19953219 -0.09465989
# sample estimates:
# mean of x mean of y 
# 0.5275780 0.6746741 
wilcox.test(top3_expression_level_vlnPlot_data[, "het_imOSN"], top3_expression_level_vlnPlot_data[, "homo_imOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top3_expression_level_vlnPlot_data[, "het_imOSN"] and top3_expression_level_vlnPlot_data[, "homo_imOSN"]
# W = 385404, p-value = 7.232e-13
# alternative hypothesis: true location shift is not equal to 0

top3_expression_level_vlnPlot_data_longer <- as.data.frame(top3_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top3_expression_level_vlnPlot_data_longer$sampletype <- factor(top3_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top3_expression_level_vlnPlot_data_longer[top3_expression_level_vlnPlot_data_longer == 0] <- NA

#top4
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
multiple_top4_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(4)(x)])
multiple_top4_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(4)(x)])
top4_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(4)(x)])
top4_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(4)(x)])
top4_list <- list(multiple_top4_expression_level_vlnPlot_data_het_mOSN, multiple_top4_expression_level_vlnPlot_data_homo_mOSN, top4_expression_level_vlnPlot_data_het_imOSN, top4_expression_level_vlnPlot_data_homo_imOSN)
top4_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top4_list, unlist), "length <-", max(lengths(top4_list))))
colnames(top4_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

t.test(top4_expression_level_vlnPlot_data[, "het_mOSN"], top4_expression_level_vlnPlot_data[, "homo_mOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top4_expression_level_vlnPlot_data[, "het_mOSN"] and top4_expression_level_vlnPlot_data[, "homo_mOSN"]
# t = -51.879, df = 9195, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -1.466952 -1.360133
# sample estimates:
# mean of x mean of y 
# 0.01251499 1.42605742
wilcox.test(top4_expression_level_vlnPlot_data[, "het_mOSN"], top4_expression_level_vlnPlot_data[, "homo_mOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top4_expression_level_vlnPlot_data[, "het_mOSN"] and top4_expression_level_vlnPlot_data[, "homo_mOSN"]
# W = 748190, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
t.test(top4_expression_level_vlnPlot_data[, "het_imOSN"], top4_expression_level_vlnPlot_data[, "homo_imOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top4_expression_level_vlnPlot_data[, "het_imOSN"] and top4_expression_level_vlnPlot_data[, "homo_imOSN"]
# t = -5.7289, df = 2966, p-value = 1.112e-08
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -0.18739718 -0.09182959
# sample estimates:
# mean of x mean of y 
# 0.3147579 0.4543713
wilcox.test(top4_expression_level_vlnPlot_data[, "het_imOSN"], top4_expression_level_vlnPlot_data[, "homo_imOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top4_expression_level_vlnPlot_data[, "het_imOSN"] and top4_expression_level_vlnPlot_data[, "homo_imOSN"]
# W = 404968, p-value = 7.111e-10
# alternative hypothesis: true location shift is not equal to 0

top4_expression_level_vlnPlot_data_longer <- as.data.frame(top4_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top4_expression_level_vlnPlot_data_longer$sampletype <- factor(top4_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top4_expression_level_vlnPlot_data_longer[top4_expression_level_vlnPlot_data_longer == 0] <- NA

#top5
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
multiple_top5_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(5)(x)])
multiple_top5_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered_multiple, 2, function(x)x[maxn(5)(x)])
top5_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(5)(x)])
top5_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(5)(x)])
top5_list <- list(multiple_top5_expression_level_vlnPlot_data_het_mOSN, multiple_top5_expression_level_vlnPlot_data_homo_mOSN, top5_expression_level_vlnPlot_data_het_imOSN, top5_expression_level_vlnPlot_data_homo_imOSN)
top5_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top5_list, unlist), "length <-", max(lengths(top5_list))))
colnames(top5_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")
t.test(top5_expression_level_vlnPlot_data[, "het_mOSN"], top5_expression_level_vlnPlot_data[, "homo_mOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top5_expression_level_vlnPlot_data[, "het_mOSN"] and top5_expression_level_vlnPlot_data[, "homo_mOSN"]
# t = -42.868, df = 9195, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -1.184252 -1.080684
# sample estimates:
# mean of x mean of y 
# 0.004006001 1.136473873
wilcox.test(top5_expression_level_vlnPlot_data[, "het_mOSN"], top5_expression_level_vlnPlot_data[, "homo_mOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top5_expression_level_vlnPlot_data[, "het_mOSN"] and top5_expression_level_vlnPlot_data[, "homo_mOSN"]
# W = 1073715, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
t.test(top5_expression_level_vlnPlot_data[, "het_imOSN"], top5_expression_level_vlnPlot_data[, "homo_imOSN"], paired = FALSE, var.equal = TRUE) # t'检验：将var.equal改为Fhet_mOSN_normalized_expression_filtered_single <- read.csv("single_functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
# Two Sample t-test
# 
# data: top5_expression_level_vlnPlot_data[, "het_imOSN"] and top5_expression_level_vlnPlot_data[, "homo_imOSN"]
# t = -4.978, df = 2966, p-value = 6.791e-07
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -0.15086911 -0.06560395
# sample estimates:
# mean of x mean of y 
# 0.2128102 0.3210467 
wilcox.test(top5_expression_level_vlnPlot_data[, "het_imOSN"], top5_expression_level_vlnPlot_data[, "homo_imOSN"], paired = F) # 两独立样本的非参数检验
# Wilcoxon rank sum test with continuity correction
# 
# data: top5_expression_level_vlnPlot_data[, "het_imOSN"] and top5_expression_level_vlnPlot_data[, "homo_imOSN"]
# W = 426282, p-value = 3.743e-07
# alternative hypothesis: true location shift is not equal to 0

top5_expression_level_vlnPlot_data_longer <- as.data.frame(top5_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top5_expression_level_vlnPlot_data_longer$sampletype <- factor(top5_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top5_expression_level_vlnPlot_data_longer[top5_expression_level_vlnPlot_data_longer == 0] <- NA

##plot

pdf("VlnPlot_top5_functional_expression_level_2.pdf", width = 7, height = 10)

p1 <- ggplot(top1_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 1")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p2 <- ggplot(top2_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 2")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p3 <- ggplot(top3_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 3")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p4 <- ggplot(top4_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 4")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p5 <- ggplot(top5_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 5")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

grid.arrange(p1, p2, p3, p4, p5, nrow = 5, ncol = 1)
dev.off()

##### plot 12.3: top5 receptors expression violin plot (all) (only_functional_receptors) ####
het_mOSN_normalized_expression_filtered <- read.csv("files_percent.mt_10_resolution_0.9/functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered <- read.csv("files_percent.mt_10_resolution_0.9/functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
het_imOSN_normalized_expression <- read.csv("files_percent.mt_10_resolution_0.9/het_imOSN_functional_normalized_expression.csv", row.names = 1, header = T)
homo_imOSN_normalized_expression <- read.csv("files_percent.mt_10_resolution_0.9/homo_imOSN_functional_normalized_expression.csv", row.names = 1, header = T)

#top1
maxn <- function(n) function(x) order(x, decreasing = TRUE)[n]
top1_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(1)(x)])
top1_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(1)(x)])
top1_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(1)(x)])
top1_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(1)(x)])
top1_list <- list(top1_expression_level_vlnPlot_data_het_mOSN, top1_expression_level_vlnPlot_data_homo_mOSN, top1_expression_level_vlnPlot_data_het_imOSN, top1_expression_level_vlnPlot_data_homo_imOSN)
top1_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top1_list, unlist), "length <-", max(lengths(top1_list))))
# do.call with the cbind function combines these vectors side-by-side into a matrix or dataframe. 
# lapply(..., "length <-", max(lengths(top1_list))): After unlisting, this applies the "length <-" function to each vector resulting from the first step. 
# lapply(top1_list, unlist): This applies the unlist function to each element in top1_list. The purpose is to convert each list within top1_list into a vector. 
colnames(top1_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

top1_expression_level_vlnPlot_data_longer <- as.data.frame(top1_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = colnames(top1_expression_level_vlnPlot_data), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top1_expression_level_vlnPlot_data_longer$sampletype <- factor(top1_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top1_expression_level_vlnPlot_data_longer[top1_expression_level_vlnPlot_data_longer == 0] <- NA

#top2
top2_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(2)(x)])
top2_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(2)(x)])
top2_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(2)(x)])
top2_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(2)(x)])
top2_list <- list(top2_expression_level_vlnPlot_data_het_mOSN, top2_expression_level_vlnPlot_data_homo_mOSN, top2_expression_level_vlnPlot_data_het_imOSN, top2_expression_level_vlnPlot_data_homo_imOSN)
top2_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top2_list, unlist), "length <-", max(lengths(top2_list))))
# do.call with the cbind function combines these vectors side-by-side into a matrix or dataframe. 
# lapply(..., "length <-", max(lengths(top2_list))): After unlisting, this applies the "length <-" function to each vector resulting from the first step. 
# lapply(top2_list, unlist): This applies the unlist function to each element in top2_list. The purpose is to convert each list within top2_list into a vector. 
colnames(top2_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

top2_expression_level_vlnPlot_data_longer <- as.data.frame(top2_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = colnames(top2_expression_level_vlnPlot_data), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top2_expression_level_vlnPlot_data_longer$sampletype <- factor(top2_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top2_expression_level_vlnPlot_data_longer[top2_expression_level_vlnPlot_data_longer == 0] <- NA

#top3
top3_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(3)(x)])
top3_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(3)(x)])
top3_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(3)(x)])
top3_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(3)(x)])
top3_list <- list(top3_expression_level_vlnPlot_data_het_mOSN, top3_expression_level_vlnPlot_data_homo_mOSN, top3_expression_level_vlnPlot_data_het_imOSN, top3_expression_level_vlnPlot_data_homo_imOSN)
top3_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top3_list, unlist), "length <-", max(lengths(top3_list))))
# do.call with the cbind function combines these vectors side-by-side into a matrix or dataframe. 
# lapply(..., "length <-", max(lengths(top3_list))): After unlisting, this applies the "length <-" function to each vector resulting from the first step. 
# lapply(top3_list, unlist): This applies the unlist function to each element in top3_list. The purpose is to convert each list within top3_list into a vector. 
colnames(top3_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

top3_expression_level_vlnPlot_data_longer <- as.data.frame(top3_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = colnames(top3_expression_level_vlnPlot_data), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top3_expression_level_vlnPlot_data_longer$sampletype <- factor(top3_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top3_expression_level_vlnPlot_data_longer[top3_expression_level_vlnPlot_data_longer == 0] <- NA

#top4
top4_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(4)(x)])
top4_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(4)(x)])
top4_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(4)(x)])
top4_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(4)(x)])
top4_list <- list(top4_expression_level_vlnPlot_data_het_mOSN, top4_expression_level_vlnPlot_data_homo_mOSN, top4_expression_level_vlnPlot_data_het_imOSN, top4_expression_level_vlnPlot_data_homo_imOSN)
top4_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top4_list, unlist), "length <-", max(lengths(top4_list))))
# do.call with the cbind function combines these vectors side-by-side into a matrix or dataframe. 
# lapply(..., "length <-", max(lengths(top4_list))): After unlisting, this applies the "length <-" function to each vector resulting from the first step. 
# lapply(top4_list, unlist): This applies the unlist function to each element in top4_list. The purpose is to convert each list within top4_list into a vector. 
colnames(top4_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

top4_expression_level_vlnPlot_data_longer <- as.data.frame(top4_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = colnames(top4_expression_level_vlnPlot_data), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top4_expression_level_vlnPlot_data_longer$sampletype <- factor(top4_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top4_expression_level_vlnPlot_data_longer[top4_expression_level_vlnPlot_data_longer == 0] <- NA

#top5
top5_expression_level_vlnPlot_data_het_mOSN <- apply(het_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(5)(x)])
top5_expression_level_vlnPlot_data_homo_mOSN <- apply(homo_mOSN_normalized_expression_filtered, 2, function(x)x[maxn(5)(x)])
top5_expression_level_vlnPlot_data_het_imOSN <- apply(het_imOSN_normalized_expression, 2, function(x)x[maxn(5)(x)])
top5_expression_level_vlnPlot_data_homo_imOSN <- apply(homo_imOSN_normalized_expression, 2, function(x)x[maxn(5)(x)])
top5_list <- list(top5_expression_level_vlnPlot_data_het_mOSN, top5_expression_level_vlnPlot_data_homo_mOSN, top5_expression_level_vlnPlot_data_het_imOSN, top5_expression_level_vlnPlot_data_homo_imOSN)
top5_expression_level_vlnPlot_data <- do.call(cbind, lapply(lapply(top5_list, unlist), "length <-", max(lengths(top5_list))))
# do.call with the cbind function combines these vectors side-by-side into a matrix or dataframe. 
# lapply(..., "length <-", max(lengths(top5_list))): After unlisting, this applies the "length <-" function to each vector resulting from the first step. 
# lapply(top5_list, unlist): This applies the unlist function to each element in top5_list. The purpose is to convert each list within top5_list into a vector. 
colnames(top5_expression_level_vlnPlot_data) <- c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN")

top5_expression_level_vlnPlot_data_longer <- as.data.frame(top5_expression_level_vlnPlot_data) %>% 
  pivot_longer( cols = colnames(top5_expression_level_vlnPlot_data), 
                names_to = 'sampletype', 
                values_to = 'expressionLevel')

top5_expression_level_vlnPlot_data_longer$sampletype <- factor(top5_expression_level_vlnPlot_data_longer$sampletype, levels = c("het_mOSN", "homo_mOSN", "het_imOSN", "homo_imOSN"))
top5_expression_level_vlnPlot_data_longer[top5_expression_level_vlnPlot_data_longer == 0] <- NA

pdf("plots_percent.mt_10_resolution_0.9/VlnPlot_top5_functional_expression_level_2_all.pdf", width = 7, height = 10)

p1 <- ggplot(top1_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 1")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p2 <- ggplot(top2_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 2")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p3 <- ggplot(top3_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 3")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p4 <- ggplot(top4_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 4")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank()) +
  xlab(NULL) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

p5 <- ggplot(top5_expression_level_vlnPlot_data_longer, aes(x = sampletype, y = expressionLevel, fill = sampletype, na.rm = T)) + 
  geom_jitter(size = 0.1, na.rm = T)+
  geom_violin(trim = T, na.rm = T)+
  scale_fill_manual(values = c("springgreen3", "lightcoral", "springgreen1", "lightpink"))+
  ggtitle("Top 5")+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_blank())+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 6)) 

grid.arrange(p1, p2, p3, p4, p5, nrow = 5, ncol = 1)
dev.off()


#### Step 13: coexpression circos plot (homo_mOSN) ####
##### plot 13.1: circos plot with coexpression pairs (all) ####
load('rdata/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"
homo_mOSN <- subset(x = trim66_harmony, idents = "mOSN_homo")
homo_mOSN_OR_expression <- homo_mOSN[grep(pattern = "^Taar|^Olfr", rownames(homo_mOSN), value = T)] 
homo_mOSN_normalized_expression <- GetAssayData(object = homo_mOSN_OR_expression, slot = "data") #(eg. “counts”, “data”, or “scale.data”). 
homo_mOSN_normalized_expression <- as.matrix(homo_mOSN_normalized_expression)
homo_mOSN_UMI <- GetAssayData(object = homo_mOSN_OR_expression, slot = "counts") #(eg. “counts”, “data”, or “scale.data”). 
homo_mOSN_UMI <- as.matrix(homo_mOSN_UMI)
all(colnames(homo_mOSN_normalized_expression) == colnames(homo_mOSN_UMI))
homo_mOSN_normalized_expression_filtered <-homo_mOSN_normalized_expression
homo_mOSN_normalized_expression_filtered[homo_mOSN_UMI < 2] <- 0

cellNumber_OR_homo_mOSN <- as.data.frame(table(colSums(homo_mOSN_normalized_expression_filtered > 0)))
colnames(cellNumber_OR_homo_mOSN) <- c("number", "homo_mOSN")

## 筛选表达多个受体的细胞
homo_mOSN_normalized_expression_filtered_multiple <- homo_mOSN_normalized_expression_filtered[, which(colSums(homo_mOSN_normalized_expression_filtered> 0)>1)]
# dim [1] 1043 8301

## 只保留有表达的基因数
homo_mOSN_normalized_expression_filtered_multiple <- homo_mOSN_normalized_expression_filtered_multiple[which(rowSums(homo_mOSN_normalized_expression_filtered_multiple) > 0), ]
# dim [1] 414 8301

multiple_receptor_name <- rownames(homo_mOSN_normalized_expression_filtered_multiple)

cellCount <- list()
n = 1
for (i in 1:413) {
  for (j in c((i+1):414)) {
    cellCount[[n]] <- as.matrix(count(colSums(homo_mOSN_normalized_expression_filtered_multiple[c(i, j), ] > 0) == 2))
    rownames(cellCount[[n]]) <- paste0(multiple_receptor_name[i], "&", multiple_receptor_name[j])
    n = n+1
  }
}

multirbind <-function(dat = list(), ...){
  if(length(dat) < 2)return(as.data.frame(dat))
  rbinddat <-dat[[1]]
  dat[[1]] <-NULL
  for(i in dat){
    rbinddat <-rbind(rbinddat, i)
  }
  return(rbinddat)
}

all_cellCount <- multirbind(cellCount)
all_cellCount_order <- as.data.frame(all_cellCount[order(all_cellCount[, 1], decreasing = T), ])
all_cellCount_order$coreceptor <- rownames(all_cellCount_order)
all_cellCount_order <- separate(data = all_cellCount_order, col = coreceptor, into = c("receptor1", "receptor2"), sep = "&")
write.csv(all_cellCount_order, "OR_coexpression_pairs_cellnumber.csv")

### circlize data preparation
#1. extract olfactory receptor information from gencode 
gencode_information <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, sep = "\t", row.names = NULL, as.is = T)

# extract taar & olfr gene
OR_information_gencode <- gencode_information[grep("^Olfr|^Taar", gencode_information$gene_name), ] #1428 rows
# only save protein_coding gene
OR_information_gencode_only_functional <- OR_information_gencode[grep("protein_coding", OR_information_gencode$gene_type), ]#1155
OR_information_gencode_only_functional <- OR_information_gencode_only_functional[-grep(pattern = "ps", OR_information_gencode_only_functional$gene_name), ]
OR_information_gencode_only_functional <- OR_information_gencode_only_functional[grep(pattern = "^chr", OR_information_gencode_only_functional$chr_name), ]

# 2. add genomic_start_relative
OR_information_gencode_only_functional$genomic_start_relative <- 1

for (i in unique(OR_information_gencode_only_functional$chr_name)) {
  OR_information_gencode_only_functional[which(OR_information_gencode_only_functional$chr_name == i), ]$genomic_start_relative <- 1:sum(OR_information_gencode_only_functional$chr_name == i)
}

# 3. add gemomic_end_relative
OR_information_gencode_only_functional$genomic_end_relative <- OR_information_gencode_only_functional$genomic_start_relative + 0.8

# 4.extract class information from lq All_OR_information
class_information <- read.csv("/lustre/home/acct-medlqian/medlqian-loop3/database/olfactory_receptor_information/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv", header = T, row.names = 1)
OR_information_gencode_only_functional$OR_class <- class_information$OR_Class[match(OR_information_gencode_only_functional$gene_name, rownames(class_information))]
OR_information_gencode_only_functional$OR_class[grep("^Taar", OR_information_gencode_only_functional$gene_name)] <- "Taar"

# 5.export necessary column
OR_information_circos_data <- OR_information_gencode_only_functional[, c("chr_name", "genomic_start_relative", "genomic_end_relative", "OR_class", "gene_name")]

# color_brewer
# http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
# https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21) 
col_for_chr_step1.1 <- col_for_chr_universal
# obtain the colors for chr used in step1.1
chr_all <- paste0("chr", c(1:19, "X", "Y"))
col_for_chr_step1.1 <- col_for_chr_step1.1[match(chr_all, unique(OR_information_circos_data$chr_name))]
col_for_chr_step1.1 <- col_for_chr_step1.1[!is.na(col_for_chr_step1.1)]

### circos plot with top 100 pairs
pdf("circosPlot_coexpression_receptor.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(OR_information_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels

# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step1.1, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the track of the different class
circos.genomicTrack(OR_information_circos_data, 
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2, 
                    bg.border = NA, 
                    track.height = 0.05, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0, 1), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop = 1, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                                 Class_I = "springgreen3", 
                                                                 Class_II = "dodgerblue3", 
                                                                 Taar = "indianred2", 
                                                                 DEFAULT = "white"), 
                                         col = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                              Class_I = "springgreen3", 
                                                              Class_II = "dodgerblue3", 
                                                              Taar = "indianred2", 
                                                              DEFAULT = "white")
                      ) #use [[]] instead of [] to obtain vector other than dataframe.
                    })

# 3. add receptor name
bed1 <- all_cellCount_functional[, 1:2]
bed1$chr <- OR_information_circos_data$chr_name[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1$start <- OR_information_circos_data$genomic_start_relative[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1$end <- OR_information_circos_data$genomic_end_relative[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1 <- bed1[, c(3, 4, 5, 1, 2)]
colnames(bed1) <- c(colnames(bed1)[1:3], "cellCount", "receptor")

bed2 <- all_cellCount_functional[, c(1, 3)]
bed2$chr <- OR_information_circos_data$chr_name[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2$start <- OR_information_circos_data$genomic_start_relative[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2$end <- OR_information_circos_data$genomic_end_relative[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2 <- bed2[, c(3, 4, 5, 1, 2)]
colnames(bed2) <- c(colnames(bed1)[1:3], "cellCount", "receptor")

bed <- rbind(bed1[1:100, ], bed2[1:100, ])
bed <- bed[!duplicated(bed$receptor), ]

circos.genomicLabels(bed, labels.column = 4, side = "outside")

# 4.add linker (coexpressed receptor)
for (i in 1:100) {
  circos.genomicLink(bed1[i, ], bed2[i, ], lwd = ifelse(bed1$cellCount[i]/1000>0.01, bed1$cellCount[i]/1000, 0.01))
}

circos.clear()
dev.off()

### circos plot with top 100 pairs only with in situ or label
pdf("circosPlot_coexpression_receptor_with_in_situ_or_label.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(OR_information_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels

# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step1.1, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the track of the different class
circos.genomicTrack(OR_information_circos_data, 
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2, 
                    bg.border = NA, 
                    track.height = 0.05, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0, 1), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop = 1, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                                 Class_I = "springgreen3", 
                                                                 Class_II = "dodgerblue3", 
                                                                 Taar = "indianred2", 
                                                                 DEFAULT = "white"), 
                                         col = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                              Class_I = "springgreen3", 
                                                              Class_II = "dodgerblue3", 
                                                              Taar = "indianred2", 
                                                              DEFAULT = "white")
                      ) #use [[]] instead of [] to obtain vector other than dataframe.
                    })

# 3. add receptor name
OR_label <- c("Olfr332", "Olfr446", "Olfr525", "Olfr874", "Olfr3", "Olfr549", "Olfr959", "Taar2", "Taar6")

OR_label_circos_data <- read.csv("OR_coexpression_pairs_cellnumber.csv", header = T) %>%
  separate(col = X, into = c("receptor1", "receptor2"), sep = "&") %>% 
  filter(receptor1 %in% OR_label) %>% 
  filter(receptor2 %in% OR_label)

bed1_OR_label <- OR_label_circos_data[, c(1, 3)]
bed1_OR_label$chr <- OR_information_circos_data$chr_name[match(bed1_OR_label$receptor1, OR_information_circos_data$gene_name)]
bed1_OR_label$start <- OR_information_circos_data$genomic_start_relative[match(bed1_OR_label$receptor1, OR_information_circos_data$gene_name)]
bed1_OR_label$end <- OR_information_circos_data$genomic_end_relative[match(bed1_OR_label$receptor1, OR_information_circos_data$gene_name)]
bed1_OR_label <- bed1_OR_label[, c(3, 4, 5, 1, 2)]
colnames(bed1_OR_label) <- c(colnames(bed1_OR_label)[1:3], "receptor", "cellCount")

bed2_OR_label <- OR_label_circos_data[, c(2, 3)]
bed2_OR_label$chr <- OR_information_circos_data$chr_name[match(bed2_OR_label$receptor2, OR_information_circos_data$gene_name)]
bed2_OR_label$start <- OR_information_circos_data$genomic_start_relative[match(bed2_OR_label$receptor2, OR_information_circos_data$gene_name)]
bed2_OR_label$end <- OR_information_circos_data$genomic_end_relative[match(bed2_OR_label$receptor2, OR_information_circos_data$gene_name)]
bed2_OR_label <- bed2_OR_label[, c(3, 4, 5, 1, 2)]
colnames(bed2_OR_label) <- c(colnames(bed1_OR_label)[1:3], "receptor", "cellCount")

bed_OR_label <- rbind(bed1_OR_label, bed2_OR_label)
bed_OR_label <- bed_OR_label[!duplicated(bed_OR_label$receptor), ] %>% na.omit()


circos.genomicLabels(bed_OR_label, labels.column = 4, side = "outside")

# 4.add linker (coexpressed receptor)

bed1 <- all_cellCount_functional[, 1:2]
bed1$chr <- OR_information_circos_data$chr_name[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1$start <- OR_information_circos_data$genomic_start_relative[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1$end <- OR_information_circos_data$genomic_end_relative[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1 <- bed1[, c(3, 4, 5, 1, 2)]
colnames(bed1) <- c(colnames(bed1)[1:3], "cellCount", "receptor")

bed2 <- all_cellCount_functional[, c(1, 3)]
bed2$chr <- OR_information_circos_data$chr_name[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2$start <- OR_information_circos_data$genomic_start_relative[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2$end <- OR_information_circos_data$genomic_end_relative[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2 <- bed2[, c(3, 4, 5, 1, 2)]
colnames(bed2) <- c(colnames(bed1)[1:3], "cellCount", "receptor")

bed <- rbind(bed1[1:100, ], bed2[1:100, ])
bed <- bed[!duplicated(bed$receptor), ]

for (i in 1:100) {
  circos.genomicLink(bed1[i, ], bed2[i, ], lwd = ifelse(bed1$cellCount[i]/1000>0.01, bed1$cellCount[i]/1000, 0.01))
}

circos.clear()
dev.off()

### only verified ors
pdf("circosPlot_coexpression_receptor_only_in_situ_ors.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(OR_information_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels

# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step1.1, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the track of the different class
circos.genomicTrack(OR_information_circos_data, 
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2, 
                    bg.border = NA, 
                    track.height = 0.05, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0, 1), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop = 1, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                                 Class_I = "springgreen3", 
                                                                 Class_II = "dodgerblue3", 
                                                                 Taar = "indianred2", 
                                                                 DEFAULT = "white"), 
                                         col = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                              Class_I = "springgreen3", 
                                                              Class_II = "dodgerblue3", 
                                                              Taar = "indianred2", 
                                                              DEFAULT = "white")
                      ) #use [[]] instead of [] to obtain vector other than dataframe.
                    })

# 3. add receptor name
OR_label <- c("Olfr332", "Olfr446", "Olfr525", "Olfr874", "Olfr3", "Olfr549", "Olfr959", "Taar2", "Taar6")

OR_label_circos_data <- read.csv("OR_coexpression_pairs_cellnumber.csv", header = T) %>%
  separate(col = X, into = c("receptor1", "receptor2"), sep = "&") %>% 
  filter(receptor1 %in% OR_label) %>% 
  filter(receptor2 %in% OR_label)


circos.genomicLabels(bed_OR_label, labels.column = 4, side = "outside")

# 4.add linker (coexpressed receptor)

for (i in 1:nrow(bed1_OR_label)) {
  circos.genomicLink(bed1_OR_label[i, ], bed2_OR_label[i, ], lwd = ifelse(bed1_OR_label$cellCount[i]/1000>0.01, bed1_OR_label$cellCount[i]/1000, 0.01))
}

circos.clear()
dev.off()

##### plot 13.2: circos plot of Taar and Olfr coexpression pairs ####
# 判断Taar跟OR是否共标，筛选出共标细胞
homo_mOSN_normalized_expression_filtered <- read.csv(file = "files_percent.mt_10_resolution_0.9/multiple_functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = 1)

homo_mOSN_normalized_expression_filtered <- homo_mOSN_normalized_expression_filtered [which(rowSums(homo_mOSN_normalized_expression_filtered) > 0), ]
# dim [1] 412 8297

Taar_Olfr_coexpression <- homo_mOSN_normalized_expression_filtered[, apply(homo_mOSN_normalized_expression_filtered[grep("^Taar", rownames(homo_mOSN_normalized_expression_filtered)), ], 2, sum)>0 & 
                                                                     apply(homo_mOSN_normalized_expression_filtered[grep("^Olfr", rownames(homo_mOSN_normalized_expression_filtered)), ], 2, sum)>0]
Taar_expression <- Taar_Olfr_coexpression[grep("^Taar", rownames(Taar_Olfr_coexpression)), ]
Olfr_expression <- Taar_Olfr_coexpression[grep("^Olfr", rownames(Taar_Olfr_coexpression)), ]

Taar_expression_name <- rownames(Taar_expression)
Olfr_expression_name <- rownames(Olfr_expression)

cellCount <- list()
a = 0
b = 1

for (i in 1:6) {
  for (j in 1:406) {
    for (n in 1:266) {
      if (Taar_expression[i, n]>0 & Olfr_expression[j, n]>0) {
        a <- a+1
      } 
    }
    cellCount[[b]] <- a
    names(cellCount[[b]]) <- paste0(Taar_expression_name[i], "&", Olfr_expression_name[j])
    a = 0
    b = b+1
  }
}

all_cellCount <- as.data.frame(multivector(cellCount))
all_cellCount_order <- arrange(all_cellCount, desc(`multivector(cellCount)`))
colnames(all_cellCount_order) <- c("cellnumber")
all_cellCount_order$coreceptor <- rownames(all_cellCount_order)
all_cellCount_order <- separate(data = all_cellCount_order, col = coreceptor, into = c("receptor1", "receptor2"), sep = "&")
write.csv(all_cellCount_order, "file/Taar_and_olfr_coexpression_cell_number_homo_multiple_mOSN.csv")
all_cellCount_order <- read.csv(file = "files_percent.mt_10_resolution_0.9/Taar_and_olfr_coexpression_cell_number_homo_multiple_mOSN.csv", header = T, row.names = 1)

### circlize data preparation
#1. extract olfactory receptor information from gencode 
gencode_information <- read.delim("/lustre/home/acct-medlqian/share/Database/Reference/Mus_musculus/GRCm38_mm10/Annotation_gtf_files/Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt", header = T, sep = "\t", row.names = NULL, as.is = T)

# extract taar & olfr gene
OR_information_gencode <- gencode_information[grep("^Olfr|^Taar", gencode_information$gene_name), ] #1428 rows
# only save protein_coding gene
OR_information_gencode_only_functional <- OR_information_gencode[grep("protein_coding", OR_information_gencode$gene_type), ]#1155
OR_information_gencode_only_functional <- OR_information_gencode_only_functional[-grep(pattern = "ps", OR_information_gencode_only_functional$gene_name), ]
OR_information_gencode_only_functional <- OR_information_gencode_only_functional[grep(pattern = "^chr", OR_information_gencode_only_functional$chr_name), ]

# 2. add genomic_start_relative
OR_information_gencode_only_functional$genomic_start_relative <- 1

for (i in unique(OR_information_gencode_only_functional$chr_name)) {
  OR_information_gencode_only_functional[which(OR_information_gencode_only_functional$chr_name == i), ]$genomic_start_relative <- 1:sum(OR_information_gencode_only_functional$chr_name == i)
}

# 3. add gemomic_end_relative
OR_information_gencode_only_functional$genomic_end_relative <- OR_information_gencode_only_functional$genomic_start_relative + 0.8

# 4.extract class information from lq All_OR_information
class_information <- read.csv("/lustre/home/acct-medlqian/medlqian-loop3/database/olfactory_receptor_information/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25_New.csv", header = T, row.names = 1)
OR_information_gencode_only_functional$OR_class <- class_information$OR_Class[match(OR_information_gencode_only_functional$gene_name, rownames(class_information))]
OR_information_gencode_only_functional$OR_class[grep("^Taar", OR_information_gencode_only_functional$gene_name)] <- "Taar"

# 5.export necessary column
OR_information_circos_data <- OR_information_gencode_only_functional[, c("chr_name", "genomic_start_relative", "genomic_end_relative", "OR_class", "gene_name")]

# color_brewer
# http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html
# https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21) 
col_for_chr_step1.1 <- col_for_chr_universal
# obtain the colors for chr used in step1.1
chr_all <- paste0("chr", c(1:19, "X", "Y"))
col_for_chr_step1.1 <- col_for_chr_step1.1[match(chr_all, unique(OR_information_circos_data$chr_name))]
col_for_chr_step1.1 <- col_for_chr_step1.1[!is.na(col_for_chr_step1.1)]

### circos plot with top 100 pairs
pdf("circosPlot_coexpression_taar_and_or_receptor_pairs.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(OR_information_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels

# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step1.1, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the track of the different class
circos.genomicTrack(OR_information_circos_data, 
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2, 
                    bg.border = NA, 
                    track.height = 0.05, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0, 1), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop = 1, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                                 Class_I = "springgreen3", 
                                                                 Class_II = "dodgerblue3", 
                                                                 Taar = "indianred2", 
                                                                 DEFAULT = "white"), 
                                         col = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                              Class_I = "springgreen3", 
                                                              Class_II = "dodgerblue3", 
                                                              Taar = "indianred2", 
                                                              DEFAULT = "white")
                      ) #use [[]] instead of [] to obtain vector other than dataframe.
                    })

# 3. add receptor name
bed1 <- all_cellCount_order[, 1:2]
bed1$chr <- OR_information_circos_data$chr_name[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1$start <- OR_information_circos_data$genomic_start_relative[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1$end <- OR_information_circos_data$genomic_end_relative[match(bed1$receptor1, OR_information_circos_data$gene_name)]
bed1 <- bed1[, c(3, 4, 5, 1, 2)]
colnames(bed1) <- c(colnames(bed1)[1:3], "cellCount", "receptor")

bed2 <- all_cellCount_order[, c(1, 3)]
bed2$chr <- OR_information_circos_data$chr_name[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2$start <- OR_information_circos_data$genomic_start_relative[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2$end <- OR_information_circos_data$genomic_end_relative[match(bed2$receptor2, OR_information_circos_data$gene_name)]
bed2 <- bed2[, c(3, 4, 5, 1, 2)]
colnames(bed2) <- c(colnames(bed1)[1:3], "cellCount", "receptor")

bed <- rbind(bed1[1:100, ], bed2[1:100, ])
bed <- bed[!duplicated(bed$receptor), ]

circos.genomicLabels(bed, labels.column = 5, side = "outside")

# 4.add linker (coexpressed receptor)
for (i in 1:100) {
  circos.genomicLink(bed1[i, ], bed2[i, ], lwd = ifelse(bed1$cellCount[i]/1000>0.01, bed1$cellCount[i]/1000, 0.01))
}

circos.clear()
dev.off()

### circos plot with top 100 pairs only with in situ or label
pdf("circosPlot_coexpression_taar_and_or_pairs_with_in_situ_or_label.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(OR_information_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels

# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step1.1, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the track of the different class
circos.genomicTrack(OR_information_circos_data, 
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2, 
                    bg.border = NA, 
                    track.height = 0.05, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0, 1), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop = 1, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                                 Class_I = "springgreen3", 
                                                                 Class_II = "dodgerblue3", 
                                                                 Taar = "indianred2", 
                                                                 DEFAULT = "white"), 
                                         col = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                              Class_I = "springgreen3", 
                                                              Class_II = "dodgerblue3", 
                                                              Taar = "indianred2", 
                                                              DEFAULT = "white")
                      ) #use [[]] instead of [] to obtain vector other than dataframe.
                    })

# 3. add receptor name
OR_label <- c("Olfr332", "Olfr446", "Olfr525", "Olfr874", "Olfr3", "Olfr549", "Olfr959", "Taar2", "Taar6")

OR_label_circos_data <- read.csv("file/Taar_and_olfr_coexpression_cell_number_homo_multiple_mOSN.csv", header = T) %>%
  separate(col = X, into = c("receptor1", "receptor2"), sep = "&") %>% 
  filter(receptor1 %in% OR_label) %>% 
  filter(receptor2 %in% OR_label)

bed1_OR_label <- OR_label_circos_data[, c(1, 3)]
bed1_OR_label$chr <- OR_information_circos_data$chr_name[match(bed1_OR_label$receptor1, OR_information_circos_data$gene_name)]
bed1_OR_label$start <- OR_information_circos_data$genomic_start_relative[match(bed1_OR_label$receptor1, OR_information_circos_data$gene_name)]
bed1_OR_label$end <- OR_information_circos_data$genomic_end_relative[match(bed1_OR_label$receptor1, OR_information_circos_data$gene_name)]
bed1_OR_label <- bed1_OR_label[, c(3, 4, 5, 1, 2)]
colnames(bed1_OR_label) <- c(colnames(bed1_OR_label)[1:3], "receptor", "cellCount")

bed2_OR_label <- OR_label_circos_data[, c(2, 3)]
bed2_OR_label$chr <- OR_information_circos_data$chr_name[match(bed2_OR_label$receptor2, OR_information_circos_data$gene_name)]
bed2_OR_label$start <- OR_information_circos_data$genomic_start_relative[match(bed2_OR_label$receptor2, OR_information_circos_data$gene_name)]
bed2_OR_label$end <- OR_information_circos_data$genomic_end_relative[match(bed2_OR_label$receptor2, OR_information_circos_data$gene_name)]
bed2_OR_label <- bed2_OR_label[, c(3, 4, 5, 1, 2)]
colnames(bed2_OR_label) <- c(colnames(bed1_OR_label)[1:3], "receptor", "cellCount")

bed_OR_label <- rbind(bed1_OR_label, bed2_OR_label)
bed_OR_label <- bed_OR_label[!duplicated(bed_OR_label$receptor), ] %>% na.omit()


circos.genomicLabels(bed_OR_label, labels.column = 4, side = "outside")

# 4.add linker (coexpressed receptor)

for (i in 1:100) {
  circos.genomicLink(bed1[i, ], bed2[i, ], lwd = ifelse(bed1$cellCount[i]/1000>0.01, bed1$cellCount[i]/1000, 0.01))
}

circos.clear()
dev.off()

### circos plot with verified taar and or pairs
pdf("circosPlot_coexpression_only_verified_taar_and_or_pairs.pdf", width = 10, height = 10)
# change the basic settings of sectors and tracks
circos.par(start.degree = 90, 
           gap.after = 0.5, # gap.after to adjust the gap between different sectors
           track.margin = c(0.008, 0.008), # set the bottom and top margin between tracks
           cell.padding = c(0, 0, 0, 0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 

# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step1.1_3, plotType = "labels")
circos.initializeWithIdeogram(OR_information_circos_data, 
                              plotType = c("labels")) # plotType = c("axis", "labels") to show axis and labels

# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step1.1, 
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the track of the different class
circos.genomicTrack(OR_information_circos_data, 
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2, 
                    bg.border = NA, 
                    track.height = 0.05, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0, 1), 
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value, 
                                         ytop = 1, ybottom = 0, # ytop.column = 1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                                 Class_I = "springgreen3", 
                                                                 Class_II = "dodgerblue3", 
                                                                 Taar = "indianred2", 
                                                                 DEFAULT = "white"), 
                                         col = hutils::Switch(value[[1]], # use [[]] instead of [] to obtain vector other than dataframe.
                                                              Class_I = "springgreen3", 
                                                              Class_II = "dodgerblue3", 
                                                              Taar = "indianred2", 
                                                              DEFAULT = "white")
                      ) #use [[]] instead of [] to obtain vector other than dataframe.
                    })

# 3. add receptor name

circos.genomicLabels(bed_OR_label, labels.column = 4, side = "outside")

# 4.add linker (coexpressed receptor)

for (i in 1:100) {
  circos.genomicLink(bed1_OR_label[i, ], bed2_OR_label[i, ], lwd = ifelse(bed1_OR_label$cellCount[i]/1000>0.01, bed1_OR_label$cellCount[i]/1000, 0.01))
}

circos.clear()
dev.off()

##### plot 13.3: circos plot of TAAR and TAAR coexpression pairs ####
# 提取共表达TAAR的细胞
homo_mOSN_normalized_expression_filtered <- read.csv(file = "multiple_functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = 1)

TAAR_homo_mOSN_normalized_expression <- homo_mOSN_normalized_expression_filtered [grep("^Taar", rownames(homo_mOSN_normalized_expression_filtered)), ]

TAAR_coexpression <- TAAR_homo_mOSN_normalized_expression[which(rowSums(TAAR_homo_mOSN_normalized_expression > 0)>0), which(colSums(TAAR_homo_mOSN_normalized_expression > 0)>1)]

cellCount <- c()
multiple_receptor_name <- rownames(TAAR_coexpression)
for (i in 1:5) {
  for (j in c((i+1):6)) {
    a <- sum(colSums(TAAR_coexpression[c(i, j), ] > 0) == 2)
    names(a) <- paste0(multiple_receptor_name[i], "&", multiple_receptor_name[j])
    cellCount <- c(cellCount, a)
  }
}

cellCount <- as.data.frame(as.matrix(cellCount))
colnames(cellCount) <- "cellNumber"

cellCount <- arrange(cellCount, desc(cellNumber))
write.csv(cellCount, "Taar_and_Taar_coexpression_cell_number_homo_multiple_mOSN.csv")

#### Step 14: bulk RNA-Seq ####
##### plot 14.1.1: correlation of receptor expression between scRNA mOSNs and bulk (only functional receptors) ####
### sc data
het_mOSN_normalized_expression_filtered <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
het_mOSN_normalized_expression_average <- apply(het_mOSN_normalized_expression_filtered, 1, mean)
homo_mOSN_normalized_expression_average <- apply(homo_mOSN_normalized_expression_filtered, 1, mean)
scRNA_mOSN_or_normalized_expression_average_ratio <- data.frame(het_mOSN_normalized_expression_average, homo_mOSN_normalized_expression_average)
write.csv(scRNA_mOSN_or_normalized_expression_average_ratio, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/scRNA_mOSN_OR_normalized_expression_average.csv")
scRNA_mOSN_or_normalized_expression_average_ratio$ratio_1 <- scRNA_mOSN_or_normalized_expression_average_ratio$homo_mOSN_normalized_expression_average/scRNA_mOSN_or_normalized_expression_average_ratio$het_mOSN_normalized_expression_average
scRNA_mOSN_or_normalized_expression_average_ratio$ratio_2 <- scRNA_mOSN_or_normalized_expression_average_ratio$homo_mOSN_normalized_expression_average/(scRNA_mOSN_or_normalized_expression_average_ratio$homo_mOSN_normalized_expression_average + scRNA_mOSN_or_normalized_expression_average_ratio$het_mOSN_normalized_expression_average)

### bulk RNA-Seq
bulk_or_DEGs <- read.csv("trim66_ko_bulk_rna_seq_results/results_functional_OR_genes_Gencode_Trim66_homo_vs_het.csv", row.names = 1, header = T)
bulk_normalized_expression <- read.csv("trim66_ko_bulk_rna_seq_results/results_Gencode_DESeq2_normalized_counts.csv", row.names = 1, header = T)
bulk_normalized_expression$bulk_het <- apply(bulk_normalized_expression[, colnames(bulk_normalized_expression)[grepl("het", colnames(bulk_normalized_expression))]], 1, mean)
bulk_normalized_expression$bulk_homo <- apply(bulk_normalized_expression[, colnames(bulk_normalized_expression)[grepl("homo", colnames(bulk_normalized_expression))]], 1, mean)
bulk_normalized_expression$bulk_ratio_1 <- bulk_normalized_expression$bulk_homo/bulk_normalized_expression$bulk_het
bulk_normalized_expression$bulk_ratio_2 <- bulk_normalized_expression$bulk_homo/(bulk_normalized_expression$bulk_het+bulk_normalized_expression$bulk_homo)

### merge
scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio <- cbind(scRNA_mOSN_or_normalized_expression_average_ratio, bulk_normalized_expression[rownames(scRNA_mOSN_or_normalized_expression_average_ratio), ])
scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na <- scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio %>% 
  drop_na(ratio_2) %>% 
  drop_na(bulk_ratio_2) # drop_na remove the rows containing na in a specified column

### plot
pdf("correlation_of_functional_genes_between_scRNA_mOSNs_and_bulk.pdf", height = 7, width = 7)
ggplot(data = scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na, aes(x = ratio_2, y = bulk_ratio_2)) + 
  geom_point() + 
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  annotate("text", x = 0.3, y = 1, 
           label = paste("Pearson Correlation: ", 
                         round(cor(scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na$ratio_2, scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na$bulk_ratio_2), 4), 
                         ", p value: ", 
                         round(cor.test(scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na$ratio_2, scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na$bulk_ratio_2)$p.value, 4), sep = ""), 
           size = 3, colour = "blue", fontface = 2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top", 
        plot.title = element_text(size = rel(2)), # adjust the title size
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.5))) +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("scRNA-Seq") + 
  ylab("bulk RNA-Seq") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
dev.off() 

##### plot 14.1.2: correlation of receptor expression between scRNA mOSNs and bulk (only DE functional receptors) ####
### DE receptor genes
bulk_or_DEGs_gene <- rownames(bulk_or_DEGs)[abs(bulk_or_DEGs$log2FoldChange)> 0.585 & bulk_or_DEGs$padj < 0.05]

scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na_DEGs <- scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na[which(rownames(scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na) %in% bulk_or_DEGs_gene), ]

pdf("correlation_of_DEGs_between_scRNA_mOSNs_and_bulk.pdf", height = 7, width = 7)
ggplot(data = scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na_DEGs, aes(x = ratio_2, y = bulk_ratio_2)) + 
  geom_point() + 
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  annotate("text", x = 0.3, y = 1, 
           label = paste("Pearson Correlation: ", 
                         round(cor(scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na_DEGs$ratio_2, scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na_DEGs$bulk_ratio_2), 4), 
                         ", p value: ", 
                         round(cor.test(scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na_DEGs$ratio_2, scRNA_mOSN_and_bulk_or_normalized_expression_average_ratio_drop_na_DEGs$bulk_ratio_2)$p.value, 4), sep = ""), 
           size = 3, colour = "blue", fontface = 2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top", 
        plot.title = element_text(size = rel(2)), # adjust the title size
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.5))) +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("scRNA-Seq") + 
  ylab("bulk RNA-Seq") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
dev.off() 

##### plot 14.2.1: violin plot: regulated OR expression downregulated and upregulated OR expression (mOSN) ####
bulk_OR_results <- read.csv("trim66_ko_bulk_rna_seq_results/results_functional_OR_genes_Gencode_Trim66_homo_vs_het.csv",header = T, row.names = 1)

downregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange < -0.585 & bulk_OR_results$padj < 0.05,])
upregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange  > 0.585 & bulk_OR_results$padj < 0.05,])
nodiff_OR_results <- na.omit(bulk_OR_results[abs(bulk_OR_results$log2FoldChange)  <= 0.585 | bulk_OR_results$padj >= 0.05,])

# import het single mOSN and homo multiple mOSN information
het_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)

downregulated_OR <- intersect(rownames(downregulated_OR_results), rownames(het_mOSN_normalized_expression_filtered))
upregulated_OR <- intersect(rownames(upregulated_OR_results), rownames(het_mOSN_normalized_expression_filtered))
nodiff_OR <- intersect(rownames(nodiff_OR_results), rownames(het_mOSN_normalized_expression_filtered))

### upregualated
het_mOSN_normalized_expression_filtered_upregulated <- het_mOSN_normalized_expression_filtered[upregulated_OR,] %>% 
  na.omit()
het_mOSN_normalized_expression_filtered_upregulated <- unlist(c(het_mOSN_normalized_expression_filtered_upregulated))
het_mOSN_normalized_expression_filtered_upregulated  <- het_mOSN_normalized_expression_filtered_upregulated[het_mOSN_normalized_expression_filtered_upregulated>0]

homo_mOSN_normalized_expression_filtered_upregulated <- homo_mOSN_normalized_expression_filtered[upregulated_OR,] %>% 
  na.omit()
homo_mOSN_normalized_expression_filtered_upregulated <- unlist(c(homo_mOSN_normalized_expression_filtered_upregulated))
homo_mOSN_normalized_expression_filtered_upregulated  <- homo_mOSN_normalized_expression_filtered_upregulated[homo_mOSN_normalized_expression_filtered_upregulated>0]

### downregulated
het_mOSN_normalized_expression_filtered_downregulated <- het_mOSN_normalized_expression_filtered[downregulated_OR,] %>% 
  na.omit()
het_mOSN_normalized_expression_filtered_downregulated <- unlist(c(het_mOSN_normalized_expression_filtered_downregulated))
het_mOSN_normalized_expression_filtered_downregulated  <- het_mOSN_normalized_expression_filtered_downregulated[het_mOSN_normalized_expression_filtered_downregulated>0]

homo_mOSN_normalized_expression_filtered_downregulated <- homo_mOSN_normalized_expression_filtered[downregulated_OR,] %>% 
  na.omit()
homo_mOSN_normalized_expression_filtered_downregulated <- unlist(c(homo_mOSN_normalized_expression_filtered_downregulated))
homo_mOSN_normalized_expression_filtered_downregulated  <- homo_mOSN_normalized_expression_filtered_downregulated[homo_mOSN_normalized_expression_filtered_downregulated>0]

### nodiff
het_mOSN_normalized_expression_filtered_nodiff <- het_mOSN_normalized_expression_filtered[nodiff_OR,] %>% 
  na.omit()
het_mOSN_normalized_expression_filtered_nodiff <- unlist(c(het_mOSN_normalized_expression_filtered_nodiff))
het_mOSN_normalized_expression_filtered_nodiff  <- het_mOSN_normalized_expression_filtered_nodiff[het_mOSN_normalized_expression_filtered_nodiff>0]

homo_mOSN_normalized_expression_filtered_nodiff <- homo_mOSN_normalized_expression_filtered[nodiff_OR,] %>% 
  na.omit()
homo_mOSN_normalized_expression_filtered_nodiff <- unlist(c(homo_mOSN_normalized_expression_filtered_nodiff))
homo_mOSN_normalized_expression_filtered_nodiff  <- homo_mOSN_normalized_expression_filtered_nodiff[homo_mOSN_normalized_expression_filtered_nodiff>0]

OR_expression <- data.frame(expression_level = c(het_mOSN_normalized_expression_filtered_upregulated, homo_mOSN_normalized_expression_filtered_upregulated,
                                                 het_mOSN_normalized_expression_filtered_downregulated, homo_mOSN_normalized_expression_filtered_downregulated,
                                                 het_mOSN_normalized_expression_filtered_nodiff, homo_mOSN_normalized_expression_filtered_nodiff),
                            diff = c(rep("het_mOSN_upregualated",length(het_mOSN_normalized_expression_filtered_upregulated)),
                                     rep("homo_mOSN_upregualated",length(homo_mOSN_normalized_expression_filtered_upregulated)),
                                     rep("het_mOSN_downregulated",length(het_mOSN_normalized_expression_filtered_downregulated)),
                                     rep("homo_mOSN_downregulated",length(homo_mOSN_normalized_expression_filtered_downregulated)),
                                     rep("het_mOSN_nodiff",length(het_mOSN_normalized_expression_filtered_nodiff)),
                                     rep("homo_mOSN_nodiff",length(homo_mOSN_normalized_expression_filtered_nodiff)))
)

OR_expression <- OR_expression %>%
  separate(col = diff, into = c("genotype", "diff"), sep = "_mOSN_", extra = "merge") %>% 
  mutate(genotype = factor(genotype, levels = c("het", "homo")),
         diff = factor(diff, levels = c("upregualated", "downregulated","nodiff")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/VlnPlot_regulated_OR_expression_level_mOSN.pdf",width = 14, height =7)

ggplot(OR_expression, aes(x = diff, y = expression_level, fill = genotype, split = genotype, na.rm = T)) + 
  scale_fill_manual(values=alpha(c("springgreen3","lightcoral"), 0.8)) +
  geom_jitter(size = 0.2, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 1)) + # jitter.width 分散的宽度，dodge.width 位置
  geom_violin(position = position_dodge(1), alpha=0.5, trim = FALSE, color = "transparent", na.rm = T) +
  facet_grid(~diff, scales='free') + # Faceting allows you to split the data into subsets and create separate plots for each subset side by side. 
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) 

dev.off()
##### plot 14.2.2: violin plot: regulated OR expression downregulated and upregulated OR expression (mOSN + imOSN) ####
bulk_OR_results <- read.csv("trim66_ko_bulk_rna_seq_results/results_functional_OR_genes_Gencode_Trim66_homo_vs_het.csv",header = T, row.names = 1)

downregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange < -0.585 & bulk_OR_results$padj < 0.05,])
upregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange  > 0.585 & bulk_OR_results$padj < 0.05,])
nodiff_OR_results <- na.omit(bulk_OR_results[abs(bulk_OR_results$log2FoldChange)  <= 0.585 | bulk_OR_results$padj >= 0.05,])

# import single mOSN and imOSN data
het_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
het_imOSN_normalized_expression <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_imOSN_functional_normalized_expression.csv", row.names = 1, header = T)
homo_imOSN_normalized_expression <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_imOSN_functional_normalized_expression.csv", row.names = 1, header = T)

downregulated_OR <- intersect(rownames(downregulated_OR_results), rownames(het_mOSN_normalized_expression_filtered))
upregulated_OR <- intersect(rownames(upregulated_OR_results), rownames(het_mOSN_normalized_expression_filtered))
nodiff_OR <- intersect(rownames(nodiff_OR_results), rownames(het_mOSN_normalized_expression_filtered))

# Define a function to pad data frames with NA rows
pad_dataframe <- function(df, max_nrow) {
  n <- nrow(df)
  if (n < max_nrow) {
    padding <- matrix(NA, nrow = max_nrow - n, ncol = ncol(df))
    colnames(padding) <- colnames(df)
    df <- rbind(df, padding)
  }
  return(df)
}

het_mOSN_normalized_expression_filtered <- pad_dataframe(het_mOSN_normalized_expression_filtered, max(nrow(het_mOSN_normalized_expression_filtered), nrow(het_imOSN_normalized_expression)))
het_imOSN_normalized_expression <- pad_dataframe(het_imOSN_normalized_expression, max(nrow(het_mOSN_normalized_expression_filtered), nrow(het_imOSN_normalized_expression)))
het_OSN_normalized_expression <- cbind(het_mOSN_normalized_expression_filtered, het_imOSN_normalized_expression)

homo_mOSN_normalized_expression_filtered <- pad_dataframe(homo_mOSN_normalized_expression_filtered, max(nrow(homo_mOSN_normalized_expression_filtered), nrow(homo_imOSN_normalized_expression)))
homo_imOSN_normalized_expression <- pad_dataframe(homo_imOSN_normalized_expression, max(nrow(homo_mOSN_normalized_expression_filtered), nrow(homo_imOSN_normalized_expression)))
homo_OSN_normalized_expression <- cbind(homo_mOSN_normalized_expression_filtered, homo_imOSN_normalized_expression)

### upregualated
het_OSN_normalized_expression_upregulated <- het_OSN_normalized_expression[upregulated_OR,] %>% 
  na.omit()
het_OSN_normalized_expression_upregulated <- unlist(c(het_OSN_normalized_expression_upregulated))
het_OSN_normalized_expression_upregulated  <- het_OSN_normalized_expression_upregulated[het_OSN_normalized_expression_upregulated>0]

homo_OSN_normalized_expression_upregulated <- homo_OSN_normalized_expression[upregulated_OR,] %>% 
  na.omit()
homo_OSN_normalized_expression_upregulated <- unlist(c(homo_OSN_normalized_expression_upregulated))
homo_OSN_normalized_expression_upregulated  <- homo_OSN_normalized_expression_upregulated[homo_OSN_normalized_expression_upregulated>0]

### downregulated
het_OSN_normalized_expression_downregulated <- het_OSN_normalized_expression[downregulated_OR,] %>% 
  na.omit()
het_OSN_normalized_expression_downregulated <- unlist(c(het_OSN_normalized_expression_downregulated))
het_OSN_normalized_expression_downregulated  <- het_OSN_normalized_expression_downregulated[het_OSN_normalized_expression_downregulated>0]

homo_OSN_normalized_expression_downregulated <- homo_OSN_normalized_expression[downregulated_OR,] %>% 
  na.omit()
homo_OSN_normalized_expression_downregulated <- unlist(c(homo_OSN_normalized_expression_downregulated))
homo_OSN_normalized_expression_downregulated  <- homo_OSN_normalized_expression_downregulated[homo_OSN_normalized_expression_downregulated>0]

### nodiff
het_OSN_normalized_expression_nodiff <- het_OSN_normalized_expression[nodiff_OR,] %>% 
  na.omit()
het_OSN_normalized_expression_nodiff <- unlist(c(het_OSN_normalized_expression_nodiff))
het_OSN_normalized_expression_nodiff  <- het_OSN_normalized_expression_nodiff[het_OSN_normalized_expression_nodiff>0]

homo_OSN_normalized_expression_nodiff <- homo_OSN_normalized_expression[nodiff_OR,] %>% 
  na.omit()
homo_OSN_normalized_expression_nodiff <- unlist(c(homo_OSN_normalized_expression_nodiff))
homo_OSN_normalized_expression_nodiff  <- homo_OSN_normalized_expression_nodiff[homo_OSN_normalized_expression_nodiff>0]

OR_expression <- data.frame(expression_level = c(het_OSN_normalized_expression_upregulated, homo_OSN_normalized_expression_upregulated,
                                                 het_OSN_normalized_expression_downregulated, homo_OSN_normalized_expression_downregulated,
                                                 het_OSN_normalized_expression_nodiff, homo_OSN_normalized_expression_nodiff),
                            diff = c(rep("het_OSN_upregualated",length(het_OSN_normalized_expression_upregulated)),
                                     rep("homo_OSN_upregualated",length(homo_OSN_normalized_expression_upregulated)),
                                     rep("het_OSN_downregulated",length(het_OSN_normalized_expression_downregulated)),
                                     rep("homo_OSN_downregulated",length(homo_OSN_normalized_expression_downregulated)),
                                     rep("het_OSN_nodiff",length(het_OSN_normalized_expression_nodiff)),
                                     rep("homo_OSN_nodiff",length(homo_OSN_normalized_expression_nodiff)))
)

OR_expression <- OR_expression %>%
  separate(col = diff, into = c("genotype", "diff"), sep = "_OSN_", extra = "merge") %>% 
  mutate(genotype = factor(genotype, levels = c("het", "homo")),
         diff = factor(diff, levels = c("upregualated", "downregulated","nodiff")))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/VlnPlot_regulated_OR_expression_level_OSN.pdf",width = 14, height =7)

ggplot(OR_expression, aes(x = diff, y = expression_level, fill = genotype, split = genotype, na.rm = T)) + 
  scale_fill_manual(values=alpha(c("springgreen3","lightcoral"), 0.8)) +
  geom_jitter(size = 0.2, position = position_jitterdodge(jitter.width = 0.4, dodge.width = 1)) + # jitter.width 分散的宽度，dodge.width 位置
  geom_violin(position = position_dodge(1), alpha=0.5, trim = FALSE, color = "transparent", na.rm = T) +
  facet_grid(~diff, scales='free') + # Faceting allows you to split the data into subsets and create separate plots for each subset side by side. 
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+ 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) 

dev.off()
##### plot 14.3.1: pie plot: regulated OR expression downregulated and upregulated OR expression cellnumber (mOSN) ####
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"
het_mOSN_cells <- subset(trim66_harmony, idents = "mOSN_het") %>% 
  ncol()
homo_mOSN_cells <- subset(trim66_harmony, idents = "mOSN_homo") %>% 
  ncol()

bulk_OR_results <- read.csv("trim66_ko_bulk_rna_seq_results/results_functional_OR_genes_Gencode_Trim66_homo_vs_het.csv",header = T, row.names = 1)
downregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange < -0.585 & bulk_OR_results$padj < 0.05,])
upregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange  > 0.585 & bulk_OR_results$padj < 0.05,])
nodiff_OR_results <- na.omit(bulk_OR_results[abs(bulk_OR_results$log2FoldChange)  <= 0.585 | bulk_OR_results$padj >= 0.05,])

# import mOSN data
het_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)

downregulated_OR <- intersect(rownames(downregulated_OR_results), rownames(het_mOSN_normalized_expression_filtered))
upregulated_OR <- intersect(rownames(upregulated_OR_results), rownames(het_mOSN_normalized_expression_filtered))
nodiff_OR <- intersect(rownames(nodiff_OR_results), rownames(het_mOSN_normalized_expression_filtered))

### upregulated
het_mOSN_normalized_expression_filtered_upregulated_cells <- het_mOSN_normalized_expression_filtered[upregulated_OR,] %>% 
  na.omit() 
het_mOSN_normalized_expression_filtered_upregulated_cells <- het_mOSN_normalized_expression_filtered_upregulated_cells[,colSums(het_mOSN_normalized_expression_filtered_upregulated_cells) > 0] %>%
  ncol()
homo_mOSN_normalized_expression_filtered_upregulated_cells <- homo_mOSN_normalized_expression_filtered[upregulated_OR,] %>% 
  na.omit() 
homo_mOSN_normalized_expression_filtered_upregulated_cells <- homo_mOSN_normalized_expression_filtered_upregulated_cells[,colSums(homo_mOSN_normalized_expression_filtered_upregulated_cells) > 0] %>%
  ncol()

### downregulated
het_mOSN_normalized_expression_filtered_downregulated_cells <- het_mOSN_normalized_expression_filtered[downregulated_OR,] %>% 
  na.omit() 
het_mOSN_normalized_expression_filtered_downregulated_cells <- het_mOSN_normalized_expression_filtered_downregulated_cells[,colSums(het_mOSN_normalized_expression_filtered_downregulated_cells) > 0] %>%
  ncol()
homo_mOSN_normalized_expression_filtered_downregulated_cells <- homo_mOSN_normalized_expression_filtered[downregulated_OR,] %>% 
  na.omit() 
homo_mOSN_normalized_expression_filtered_downregulated_cells <- homo_mOSN_normalized_expression_filtered_downregulated_cells[,colSums(homo_mOSN_normalized_expression_filtered_downregulated_cells) > 0] %>%
  ncol()

### nodiff
het_mOSN_normalized_expression_filtered_nodiff_cells <- het_mOSN_normalized_expression_filtered[nodiff_OR,] %>% 
  na.omit() 
het_mOSN_normalized_expression_filtered_nodiff_cells <- het_mOSN_normalized_expression_filtered_nodiff_cells[,colSums(het_mOSN_normalized_expression_filtered_nodiff_cells) > 0] %>%
  ncol()
homo_mOSN_normalized_expression_filtered_nodiff_cells <- homo_mOSN_normalized_expression_filtered[nodiff_OR,] %>% 
  na.omit() 
homo_mOSN_normalized_expression_filtered_nodiff_cells <- homo_mOSN_normalized_expression_filtered_nodiff_cells[,colSums(homo_mOSN_normalized_expression_filtered_nodiff_cells) > 0] %>%
  ncol()

OR_expression_cells <- data.frame(cellNumber_het = c(het_mOSN_cells, het_mOSN_normalized_expression_filtered_upregulated_cells, het_mOSN_normalized_expression_filtered_downregulated_cells, het_mOSN_normalized_expression_filtered_nodiff_cells),
                                  cellNumber_homo = c(homo_mOSN_cells, homo_mOSN_normalized_expression_filtered_upregulated_cells, homo_mOSN_normalized_expression_filtered_downregulated_cells, homo_mOSN_normalized_expression_filtered_nodiff_cells),
                                  row.names = c("all", "upregulated", "downregulated", "nodiff"))
# cellNumber_het cellNumber_homo
# all                     6223            9430
# upregulated              471            8535
# downregulated           4986            2173
# nodiff                   400            4976

# Calculate percentages for each category
OR_expression_cells <- OR_expression_cells %>% 
  mutate(percent_het = OR_expression_cells$cellNumber_het / OR_expression_cells$cellNumber_het[1] * 100,
         percent_homo = OR_expression_cells$cellNumber_homo / OR_expression_cells$cellNumber_homo[1] *100,
         category = rownames(.)) %>% 
  filter(rownames(.) != "all") 
# cellNumber_het cellNumber_homo percent_het percent_homo      category
# upregulated              471            8535    7.568697     90.50901   upregulated
# downregulated           4986            2173   80.122128     23.04348 downregulated
# nodiff                   400            4976    6.427768     52.76776        nodiff
write.csv(OR_expression_cells, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/data_for_piePlot_regulated_ORs_expressing_cellNumber_mOSN.csv")

# Function to plot a single pie chart for a given category
plot_pie <- function(category, percentages, labels, my_colors) {
  pie(percentages, labels = labels, main = paste("Pie Chart for", category), col = my_colors)
}

# Loop through categories (excluding 'all') and plot pie charts
categories <- c("upregulated", "downregulated", "nodiff")

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/piePlot_regulated_ORs_expressing_cellNumber_mOSN.pdf", paper = "a4")
par(mfrow=c(3,2)) 
for (category in categories) {
  cat_OR_expression_cells <- OR_expression_cells[OR_expression_cells$category == category, ]
  # Plot for het
  plot_pie(category, c(cat_OR_expression_cells$percent_het, 100 - cat_OR_expression_cells$percent_het), c("Het", "Others"),
           c("springgreen3", "grey"))
  # Plot for homo
  plot_pie(category, c(cat_OR_expression_cells$percent_homo, 100 - cat_OR_expression_cells$percent_homo), c("Homo", "Others"),
           c("lightcoral", "grey"))
}
dev.off()
##### plot 14.3.2: pie plot: regulated OR expression downregulated and upregulated OR expression cellnumber (mOSN + imOSN) ####
load('results_percent.mt_10_integrate_resolution_9/rdata_percent.mt_10_resolution_0.9/trim66_harmony_percent.mt_10_resolution_0.9_afterAddMetadata.rdata')
Idents(trim66_harmony) <- "celltype.genotype"
het_OSN_cells <- subset(trim66_harmony, idents = c("mOSN_het", "imOSN_het")) %>% 
  ncol()
homo_OSN_cells <- subset(trim66_harmony, idents = c("mOSN_homo", "imOSN_homo")) %>% 
  ncol()

bulk_OR_results <- read.csv("trim66_ko_bulk_rna_seq_results/results_functional_OR_genes_Gencode_Trim66_homo_vs_het.csv",header = T, row.names = 1)
downregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange < -0.585 & bulk_OR_results$padj < 0.05,])
upregulated_OR_results <- na.omit(bulk_OR_results[bulk_OR_results$log2FoldChange  > 0.585 & bulk_OR_results$padj < 0.05,])
nodiff_OR_results <- na.omit(bulk_OR_results[abs(bulk_OR_results$log2FoldChange)  <= 0.585 | bulk_OR_results$padj >= 0.05,])

# import mOSN and imOSN data
het_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_het_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
homo_mOSN_normalized_expression_filtered <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/functional_or_homo_mOSN_normalized_expression_filtered_UMI_2.csv", row.names = 1, header = T)
het_imOSN_normalized_expression <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_imOSN_functional_normalized_expression.csv", row.names = 1, header = T)
homo_imOSN_normalized_expression <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_imOSN_functional_normalized_expression.csv", row.names = 1, header = T)

# Define a function to pad data frames with NA rows
pad_dataframe <- function(df, max_nrow) {
  n <- nrow(df)
  if (n < max_nrow) {
    padding <- matrix(NA, nrow = max_nrow - n, ncol = ncol(df))
    colnames(padding) <- colnames(df)
    df <- rbind(df, padding)
  }
  return(df)
}

het_mOSN_normalized_expression_filtered <- pad_dataframe(het_mOSN_normalized_expression_filtered, max(nrow(het_mOSN_normalized_expression_filtered), nrow(het_imOSN_normalized_expression)))
het_imOSN_normalized_expression <- pad_dataframe(het_imOSN_normalized_expression, max(nrow(het_mOSN_normalized_expression_filtered), nrow(het_imOSN_normalized_expression)))
het_OSN_normalized_expression <- cbind(het_mOSN_normalized_expression_filtered, het_imOSN_normalized_expression)

homo_mOSN_normalized_expression_filtered <- pad_dataframe(homo_mOSN_normalized_expression_filtered, max(nrow(homo_mOSN_normalized_expression_filtered), nrow(homo_imOSN_normalized_expression)))
homo_imOSN_normalized_expression <- pad_dataframe(homo_imOSN_normalized_expression, max(nrow(homo_mOSN_normalized_expression_filtered), nrow(homo_imOSN_normalized_expression)))
homo_OSN_normalized_expression <- cbind(homo_mOSN_normalized_expression_filtered, homo_imOSN_normalized_expression)

downregulated_OR <- intersect(rownames(downregulated_OR_results), rownames(het_OSN_normalized_expression_filtered))
upregulated_OR <- intersect(rownames(upregulated_OR_results), rownames(het_OSN_normalized_expression_filtered))
nodiff_OR <- intersect(rownames(nodiff_OR_results), rownames(het_OSN_normalized_expression_filtered))

downregulated_OR <- intersect(rownames(downregulated_OR_results), rownames(het_OSN_normalized_expression))
upregulated_OR <- intersect(rownames(upregulated_OR_results), rownames(het_OSN_normalized_expression))
nodiff_OR <- intersect(rownames(nodiff_OR_results), rownames(het_OSN_normalized_expression))

### upregulated
het_OSN_normalized_expression_upregulated_cells <- het_OSN_normalized_expression[upregulated_OR,] %>% 
  na.omit() 
het_OSN_normalized_expression_upregulated_cells <- het_OSN_normalized_expression_upregulated_cells[,colSums(het_OSN_normalized_expression_upregulated_cells) > 0] %>%
  ncol()
homo_OSN_normalized_expression_upregulated_cells <- homo_OSN_normalized_expression[upregulated_OR,] %>% 
  na.omit() 
homo_OSN_normalized_expression_upregulated_cells <- homo_OSN_normalized_expression_upregulated_cells[,colSums(homo_OSN_normalized_expression_upregulated_cells) > 0] %>%
  ncol()

### downregulated
het_OSN_normalized_expression_downregulated_cells <- het_OSN_normalized_expression[downregulated_OR,] %>% 
  na.omit() 
het_OSN_normalized_expression_downregulated_cells <- het_OSN_normalized_expression_downregulated_cells[,colSums(het_OSN_normalized_expression_downregulated_cells) > 0] %>%
  ncol()
homo_OSN_normalized_expression_downregulated_cells <- homo_OSN_normalized_expression[downregulated_OR,] %>% 
  na.omit() 
homo_OSN_normalized_expression_downregulated_cells <- homo_OSN_normalized_expression_downregulated_cells[,colSums(homo_OSN_normalized_expression_downregulated_cells) > 0] %>%
  ncol()

### nodiff
het_OSN_normalized_expression_nodiff_cells <- het_OSN_normalized_expression[nodiff_OR,] %>% 
  na.omit() 
het_OSN_normalized_expression_nodiff_cells <- het_OSN_normalized_expression_nodiff_cells[,colSums(het_OSN_normalized_expression_nodiff_cells) > 0] %>%
  ncol()
homo_OSN_normalized_expression_nodiff_cells <- homo_OSN_normalized_expression[nodiff_OR,] %>% 
  na.omit() 
homo_OSN_normalized_expression_nodiff_cells <- homo_OSN_normalized_expression_nodiff_cells[,colSums(homo_OSN_normalized_expression_nodiff_cells) > 0] %>%
  ncol()

OR_expression_cells <- data.frame(cellNumber_het = c(het_OSN_cells, het_OSN_normalized_expression_upregulated_cells, het_OSN_normalized_expression_downregulated_cells, het_OSN_normalized_expression_nodiff_cells),
                                  cellNumber_homo = c(homo_OSN_cells, homo_OSN_normalized_expression_upregulated_cells, homo_OSN_normalized_expression_downregulated_cells, homo_OSN_normalized_expression_nodiff_cells),
                                  row.names = c("all", "upregulated", "downregulated", "nodiff"))
# cellNumber_het cellNumber_homo
# all                     6729           12561
# upregulated              557           10859
# downregulated           5374            4304
# nodiff                   480            6094

# Calculate percentages for each category
OR_expression_cells <- OR_expression_cells %>% 
  mutate(percent_het = OR_expression_cells$cellNumber_het / OR_expression_cells$cellNumber_het[1] * 100,
         percent_homo = OR_expression_cells$cellNumber_homo / OR_expression_cells$cellNumber_homo[1] *100,
         category = rownames(.)) %>% 
  filter(rownames(.) != "all") 
# cellNumber_het cellNumber_homo percent_het percent_homo      category
# upregulated              557           10859    8.277604     86.45012   upregulated
# downregulated           5374            4304   79.863278     34.26479 downregulated
# nodiff                   480            6094    7.133304     48.51525        nodiff
write.csv(OR_expression_cells, "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/data_for_piePlot_regulated_ORs_expressing_cellNumber_OSN.csv")

# Function to plot a single pie chart for a given category
plot_pie <- function(category, percentages, labels, my_colors) {
  pie(percentages, labels = labels, main = paste("Pie Chart for", category), col = my_colors)
}

# Loop through categories (excluding 'all') and plot pie charts
categories <- c("upregulated", "downregulated", "nodiff")

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/piePlot_regulated_ORs_expressing_cellNumber_OSN.pdf", paper = "a4")
par(mfrow=c(3,2)) 
for (category in categories) {
  cat_OR_expression_cells <- OR_expression_cells[OR_expression_cells$category == category, ]
  # Plot for het
  plot_pie(category, c(cat_OR_expression_cells$percent_het, 100 - cat_OR_expression_cells$percent_het), c("Het", "Others"),
           c("springgreen3", "grey"))
  # Plot for homo
  plot_pie(category, c(cat_OR_expression_cells$percent_homo, 100 - cat_OR_expression_cells$percent_homo), c("Homo", "Others"),
           c("lightcoral", "grey"))
}
dev.off()

##### Plot 16.3: box plot: downregulated and upregulated OR expression cell number (mOSN) ####
het_mOSN_single_functional_Taar_OR_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_single_functional_Taar_OR_expression_level.csv", row.names = 1, header = T)
het_mOSN_multiple_functional_Taar_OR_level <- read.csv(file = "results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/het_mOSN_multiple_functional_Taar_OR_expression_level.csv", row.names = 1, header = T)
homo_mOSN_single_functional_Taar_OR_level <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_mOSN_single_functional_Taar_OR_expression_level.csv", row.names = 1, header = T)
homo_mOSN_multiple_functional_Taar_OR_level <- read.csv("results_percent.mt_10_integrate_resolution_9/files_percent.mt_10_resolution_0.9/homo_mOSN_multiple_functional_Taar_OR_expression_level.csv", row.names = 1, header = T)

upregualated_OR_het_single <- het_mOSN_single_functional_Taar_OR_level[rownames(upregulated_OR), ]
upregualated_OR_het_single <- upregualated_OR_het_single$cellNumber
upregualated_OR_het_single <- upregualated_OR_het_single[which(upregualated_OR_het_single>0)]

downregulated_OR_het_single <- het_mOSN_single_functional_Taar_OR_level[rownames(downregulated_OR), ]
downregulated_OR_het_single <- downregulated_OR_het_single$cellNumber
downregulated_OR_het_single <- downregulated_OR_het_single[which(downregulated_OR_het_single>0)]

upregualated_OR_het_multiple <- het_mOSN_multiple_functional_Taar_OR_level[rownames(upregulated_OR), ]
upregualated_OR_het_multiple <- upregualated_OR_het_multiple$cellNumber
upregualated_OR_het_multiple <- upregualated_OR_het_multiple[which(upregualated_OR_het_multiple>0)]

downregulated_OR_het_multiple <- het_mOSN_multiple_functional_Taar_OR_level[rownames(downregulated_OR), ]
downregulated_OR_het_multiple <- downregulated_OR_het_multiple$cellNumber
downregulated_OR_het_multiple <- downregulated_OR_het_multiple[which(downregulated_OR_het_multiple>0)]

upregualated_OR_homo_single <- homo_mOSN_single_functional_Taar_OR_level[rownames(upregulated_OR), ]
upregualated_OR_homo_single <- upregualated_OR_homo_single$cellNumber
upregualated_OR_homo_single <- upregualated_OR_homo_single[which(upregualated_OR_homo_single>0)]

downregulated_OR_homo_single <- homo_mOSN_single_functional_Taar_OR_level[rownames(downregulated_OR), ]
downregulated_OR_homo_single <- downregulated_OR_homo_single$cellNumber
downregulated_OR_homo_single <- downregulated_OR_homo_single[which(downregulated_OR_homo_single>0)]

upregualated_OR_homo_multiple <- homo_mOSN_multiple_functional_Taar_OR_level[rownames(upregulated_OR), ]
upregualated_OR_homo_multiple <- upregualated_OR_homo_multiple$cellNumber
upregualated_ORupregualated_OR_homo_multiple_homo <- upregualated_OR_homo_multiple[which(upregualated_OR_homo_multiple>0)]

downregulated_OR_homo_multiple <- homo_mOSN_multiple_functional_Taar_OR_level[rownames(downregulated_OR), ]
downregulated_OR_homo_multiple <- downregulated_OR_homo_multiple$cellNumber
downregulated_OR_homo_multiple <- downregulated_OR_homo_multiple[which(downregulated_OR_homo_multiple>0)]

OR_cell_number <- data.frame(cellNumber = c(upregualated_OR_het_single, downregulated_OR_het_single, 
                                            upregualated_OR_het_multiple, downregulated_OR_het_multiple, 
                                            upregualated_OR_homo_single, downregulated_OR_homo_single, 
                                            upregualated_OR_homo_multiple, downregulated_OR_homo_multiple), 
                             diff = c(rep("upregualated_OR_cell_number_het_single", length(upregualated_OR_het_single)), 
                                      rep("downregulated_OR_cell_number_het_single", length(downregulated_OR_het_single)), 
                                      rep("upregualated_OR_cell_number_het_multiple", length(upregualated_OR_het_multiple)), 
                                      rep("downregulated_OR_cell_number_het_multiple", length(downregulated_OR_het_multiple)), 
                                      rep("upregualated_OR_cell_number_homo_single", length(upregualated_OR_homo_single)), 
                                      rep("downregulated_OR_cell_number_homo_single", length(downregulated_OR_homo_single)), 
                                      rep("upregualated_OR_cell_number_homo_multiple", length(upregualated_OR_homo_multiple)), 
                                      rep("downregulated_OR_cell_number_homo_multiple", length(downregulated_OR_homo_multiple))))

OR_cell_number_separate <- separate(data = OR_cell_number, col = diff, into = c("upOrDown", "genotype"), sep = "_OR_cell_number_")

OR_cell_number <- cbind(OR_cell_number, OR_cell_number_separate[2:3])

OR_cell_number$genotype <- factor(OR_cell_number$genotype, levels = c("het_single", "het_multiple", "homo_single", "homo_multiple"))
OR_cell_number$upOrDown <- factor(OR_cell_number$upOrDown, levels = c("upregualated", "downregulated"))

pdf("results_percent.mt_10_integrate_resolution_9/plots_percent.mt_10_resolution_0.9/BoxPlot_cellNumber_downOrUp_regulated_OR_mOSN_modified.pdf", width = 10, height = 7)

down <- ggplot(OR_cell_number, aes(x = upOrDown, y = cellNumber, fill = genotype, split = genotype, na.rm = T)) +
  scale_fill_manual(values = alpha(c("springgreen3", "springgreen1", "lightcoral", "lightpink"), 0.8)) +
  geom_jitter(size = 0.2, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1))+
  geom_boxplot(position = position_dodge(1), outlier.size = 0.2, alpha = 0.5) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 30, 10)) +
  theme_classic() +
  labs(x = "", y = "") +
  theme(legend.position = "bottom") +
  coord_cartesian(ylim = c(0, 30))+
  facet_grid(~upOrDown, scales = 'free') + 
  theme(strip.text = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

upper <- ggplot(OR_cell_number, aes(x = upOrDown, y = cellNumber, fill = genotype, split = genotype, na.rm = T)) +
  scale_fill_manual(values = alpha(c("springgreen3", "springgreen1", "lightcoral", "lightpink"), 0.8)) +
  geom_jitter(size = 0.2, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1)) +
  geom_boxplot(position = position_dodge(1), outlier.size = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(x = "", y = "") +
  theme(legend.position = "none") +
  coord_cartesian(ylim = c(30, 6000)) + 
  scale_y_continuous(breaks = c(100, 2000, 4000, 6000)) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank())+
  facet_grid(~upOrDown, scales = 'free') 

# install.packages("ggpubr")
# library(ggpubr)
ggarrange(upper, 
          down, 
          heights = c(1, 1), 
          widths = c(1, 1), 
          ncol = 1, 
          nrow = 2, 
          align = "hv")

dev.off()
