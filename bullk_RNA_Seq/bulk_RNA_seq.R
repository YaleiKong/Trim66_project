
## The bulk RNA-seq data include three genotypes: Trim66 knockout, Foxg1-Cre; Trim66 flox/flox, 
#and Goofy-Cre; Trim66 flox/flox mice. The code provided below illustrates the 
#analysis workflow using the Trim66 knockout bulk RNA-seq dataset as an example. 
#Identical pipelines were used for the other two genotypes.

library(geneplotter)
library(plyr)
library(dplyr)
library(genefilter)
library(topGO)
library(stringr)
library(LSD)
library(tidyr)
library(biomaRt)
library(EDASeq)
library(DESeq2)
library(gProfileR)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(gplots)
library(reshape2)
library(ggpubr)
library(ggbreak)
library(gtools) # to use mixedsort, check http://veda.cs.uiuc.edu/TCGA_classify/RWR/DRaWR/library/gtools/html/mixedsort.html
library(circlize)
library(ComplexHeatmap) # to add legend in Circos plot
library(hutils) # to use Switch, which is the vectorized form of base::switch. like ifelse is vevtorized form of if else. Check: https://www.delftstack.com/howto/r/use-a-vectorized-if-function-with-multiple-conditions-in-r/
library(readxl)
### Set variables
rm(list=ls(all=TRUE)) #clear all the variables
setwd('/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Trim66_Cbx4_project/Trim66_KO_bulk_seq/20220126_Trim66_KO_bulk_seq/DESeq2_analysis_Gencode')
baseDir <- getwd()
dataDir <- paste(baseDir, "/data", sep="")
resultDir <- paste(baseDir, "/results", sep="")
metaDir <- paste(baseDir, "/meta", sep="")
# save the R sessionInfo to a txt file for future reference
writeLines(capture.output(sessionInfo()), "Trim66_KO_RNA_seq_20220501_sessionInfo.txt")
# set the Gencode database directory
DatabaseDir <- '/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Gencode/mouse_genes'
## Load data and metadata
# as.is=T means keep all character vectors, do not convert to factors
data <- read.delim(file.path(dataDir, "all_trimmomatic_featureCounts_with_header.txt"), header=T, sep="\t",
                   row.names=1, as.is=T)
row_names <- rownames(data)
# load metadata which is the sample categories
meta <- read.delim(file.path(metaDir, "Trim66_KO_RNA_seq_design.txt"), header=T, sep="\t",
                   row.names=1, stringsAsFactors=T)
dim(meta)
summary(meta)
## Check that sample names match in both files
# "all" function to check whether all the logic are true
which(row.names(meta) %in% colnames(data))
all(row.names(meta) == colnames(data))

# load ensembl gene name and gene description file
ensembl_gene_list <- read.delim(file.path(DatabaseDir, "Genes_extracted_from_gencode.vM25.Comprehensive.gtf.file_modified_final.txt"), header=T, sep="\t",
                   row.names=NULL, as.is=T)
# add one more column with the ensembl ID without the suffix that are the version information
ensembl_gene_list$gene_id_Ensembl_remove_version_info <- str_remove(ensembl_gene_list$gene_id_Ensembl, "\\.[0-9]")
##########################################################################
    #Note: obtain protein-coding genes and pseudogenes including unprocessed_pseudogene and polymorphic_pseudogenes that OR/TAAR pseudogenes belong. Check https://www.gencodegenes.org/pages/biotypes.html
##########################################################################
ensembl_gene_list_good <- ensembl_gene_list[grepl("protein_coding|unprocessed_pseudogene|polymorphic_pseudogene", 
                                                  ensembl_gene_list$gene_type), ] # use "|" to match multiple strings in grepl, note no space between strings
dim(ensembl_gene_list_good)
#exclude 13 mitochondria genes (chrM) and other unknown chr labels. Keep only chr1, chr2, ....chrX, and chrY
sum(grepl("chrM", ensembl_gene_list_good$chr_name))
sum(!grepl("chr", ensembl_gene_list_good$chr_name))
ensembl_gene_list_good <- ensembl_gene_list_good[grepl("^chr[1-9XY]", ensembl_gene_list_good$chr_name), ]
dim(ensembl_gene_list_good)

# match the ensembl ID of the data to the ensembl ID from ensembl_gene_list_good
good_data <- data[rownames(data) %in% ensembl_gene_list_good$gene_id_Ensembl, ]
# obtain the gene names
good_data_matched_gene_name <- good_data
# make.names with unique=T parameter as there are several genes have multiple gene ID. Here 30 duplicated names, such as Ndor1, Ptp4a1, Olfr1073-ps1, Olfr290, Vmn1r216 ...
gene_names_from_list <- ensembl_gene_list_good[match(rownames(good_data_matched_gene_name), ensembl_gene_list_good$gene_id_Ensembl), ]$gene_name
duplicated_gene_names <- make.names(unique(gene_names_from_list[duplicated(gene_names_from_list)]))
rownames(good_data_matched_gene_name) <- make.names(ensembl_gene_list_good[match(rownames(good_data_matched_gene_name), ensembl_gene_list_good$gene_id_Ensembl), ]$gene_name, unique = T)
dim(good_data_matched_gene_name)
# keep a record of matched gene ID and gene name
ensembl_genes_from_seq_data <- ensembl_gene_list_good[match(rownames(good_data), ensembl_gene_list_good$gene_id_Ensembl), ]
ensembl_genes_from_seq_data$Gene_name_in_seq_data <- make.names(ensembl_gene_list_good[match(rownames(good_data), ensembl_gene_list_good$gene_id_Ensembl), ]$gene_name, unique = T)
ensembl_genes_from_seq_data[ensembl_genes_from_seq_data$gene_name == "Ndor1", ] # for instance, gene Arf4 has two IDs.
# remove the genes with multiple gene IDs caused by make.names above (having pattern like gene.1, gene.2, gene.3 et al), as I found that those genes with .1, .2 et al show NA in the following DESeq2 results
duplicated_gene_names_pattern <- paste(duplicated_gene_names, ".[1-9]", sep = "") # matched pattern is duplicated gene name followed by "."
duplicated_gene_names_pattern_for_grepl <- paste(duplicated_gene_names_pattern, collapse = "|")
good_data_matched_gene_name <- good_data_matched_gene_name[!grepl(duplicated_gene_names_pattern_for_grepl, 
                                                                 rownames(good_data_matched_gene_name)), ] # "\\" to escape special character. check https://cran.r-project.org/web/packages/stringr/vignettes/regular-expressions.html
dim(good_data_matched_gene_name)
dim(good_data)
## Quality control (pre-normalization)

################################ PCA analysis of original counts data  ######################
# PCA (principal components analysis) is a multivariate technique that allows us to summarize 
# the systematic patterns of variations in the data. PCA takes the expresson levels for all 
# probes and transforms it in principal component space, reducing each sample into one point. 
# This allows us to separate samples according to expression variation, and identify potential 
# outliers. 

# PCA plot of first and second PCs
#prcomp to make PCA analysis, it will generate a list data type, "t" means transpose which is
#required by the function
# "x" in the prcomp is a matrix containing all the PCA parameters
pca_matrix <- prcomp(t(good_data_matched_gene_name))
names(pca_matrix)
pca_matrix <- pca_matrix$x
df_pca <- cbind(meta, pca_matrix[, c('PC1', 'PC2')])
sampletype <- meta[,'sampletype']
ggplot(df_pca, aes(PC1, PC2, color= sampletype)) +
  geom_text(aes(PC1, PC2, label=row.names(df_pca)), size=5, hjust=0.1, vjust=0.1) +
  scale_x_continuous(expand=c(0.3, 0.3))
### Sample-to-sample correlation heatmap
# Another method to identify potential outliers
# "annotation" table is needed for pheatmap, use data.frame function to generate
# it. Because it need row names
# Set color palettes for visualization tools
#"rev" function to reverse the sequence of the vector, brewer.pal function to generate a color range
heatcolors.1 <- rev(brewer.pal(6, "YlOrRd"))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "F0E442", "#0072B2",
               "#D55E00", "#CC79A7", "#000000")
annotation <- data.frame(sampletype=meta[,"sampletype"],
                         row.names=row.names(meta))
pheatmap(cor(good_data_matched_gene_name), color=heatcolors.1, cluster_rows=T, annotation=annotation,
         border_color=NA, cluster_cols=T, clustering_distance_rows="euclidean",
         clustering_distance_cols="euclidean", fontsize=10, fontsize_row=8, height=20)
#################################################################################################


#################Find  differential Expression genes #############################################
dds_data <- DESeqDataSetFromMatrix(countData=good_data_matched_gene_name, colData=meta, design= ~ sampletype)
dds_data$sampletype <- relevel(dds_data$sampletype, "het")
dds <- DESeq(dds_data)  # type vignette("DESeq2") for detailed information
#what are the size factors for each samples?
sizeFactors(dds)
# what are dispersion factors for every gene. They reflect within-group variability and are used for statistical analysis.
dispersions(dds)
# what are data stored in dds
mcols(dds, use.names = TRUE) #for example, use mcols(dds)$maxCooks to extract the data from column maxCooks
mcols(mcols(dds), use.names = T) # explain what each columns mean
# to obtain the original counts or normalized counts
head(counts(dds))
normalized_counts <- DESeq2::counts(dds, normalized=TRUE) #normalized counts This is just dividing each column of counts(dds) by sizeFactors(dds), which is equal to sweep(counts(dds), 2, sizeFactors(dds), "/") where 2 means column wise
plotCounts(dds, gene="Olfr51", intgroup=c("sampletype")) #plotCounts to inspect and plot normalized counts of genes you are interested
barplot(normalized_counts[rownames(normalized_counts)=="Olfr51", ],las=2, main=rownames(normalized_counts)[rownames(normalized_counts)=="Olfr51"]  )


#If there are more than 2 levels for the variable (as is the case
# in this analysis) results will extract the results table for a comparison of the last level over the first level.
# set independentFiltering to FALSE so that there will be less NA values of padj (see below for padj). Especially for
# lowly expressed genes such as olfactory receptors. The independent filtering is designed only to filter out low count genes to the extent that they are not enriched with small p-values.
res_all <- results(dds, independentFiltering = FALSE)
# add "gene_type", "ensembl ID", "chr", and "genomic_strand" information from above ensembl_genes_from_seq_data variable in res_all
res_all$gene_type <- ensembl_genes_from_seq_data[match(rownames(res_all), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_type
res_all$gene_id_Ensembl <- ensembl_genes_from_seq_data[match(rownames(res_all), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_id_Ensembl
res_all$gene_id_Ensembl_remove_version_info <- ensembl_genes_from_seq_data[match(rownames(res_all), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_id_Ensembl_remove_version_info
res_all$chr_name <- ensembl_genes_from_seq_data[match(rownames(res_all), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$chr_name
res_all$genomic_strand <- ensembl_genes_from_seq_data[match(rownames(res_all), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_strand
res_all$genomic_start <- ensembl_genes_from_seq_data[match(rownames(res_all), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_start
res_all$genomic_end <- ensembl_genes_from_seq_data[match(rownames(res_all), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_end
head(res_all)  # In the output, baseMean = the average of the normalized counts taken over all samples of level 1. Here it is WT as we relevel wild type as level 1.
# log2FoldChange = log2 fold change between the groups. E.g. value 2 means that the expression has increased 4-fold
# lfcSE = standard error of the log2FoldChange estimate
# stat = Wald or LRT statistic. Here LRT statistic
# pvalue = Wald or LRT test p-value. Here LRT p value
# padj = Benjamini-Hochberg adjusted p-value. DESeq2 uses the so-called Benjamini-Hochberg (BH) adjustment which is padj in res_homo_WT for example. This is more 
# suitable than pvalue

#################################################################################################################################
#### write to file
write.csv(res_all, file = paste(baseDir, "/results/results_Gencode_Trim66_homo_vs_het.csv", sep=''))
# save the raw counts
original_counts <- DESeq2::counts(dds)
original_counts <- as.data.frame(original_counts)
original_counts$gene_type <- ensembl_genes_from_seq_data[match(rownames(original_counts), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_type
original_counts$gene_id_Ensembl <- ensembl_genes_from_seq_data[match(rownames(original_counts), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_id_Ensembl
original_counts$gene_id_Ensembl_remove_version_info <- ensembl_genes_from_seq_data[match(rownames(original_counts), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_id_Ensembl_remove_version_info
original_counts$chr_name <- ensembl_genes_from_seq_data[match(rownames(original_counts), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$chr_name
original_counts$genomic_strand <- ensembl_genes_from_seq_data[match(rownames(original_counts), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_strand
original_counts$genomic_start <- ensembl_genes_from_seq_data[match(rownames(original_counts), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_start
original_counts$genomic_end <- ensembl_genes_from_seq_data[match(rownames(original_counts), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_end
write.csv(original_counts, file = paste(baseDir, "/results/results_Gencode_raw_counts.csv", sep=''))
# save normalized counts
normalized_counts_dataframe <- as.data.frame(normalized_counts)
normalized_counts_dataframe$gene_type <- ensembl_genes_from_seq_data[match(rownames(normalized_counts_dataframe), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_type
normalized_counts_dataframe$gene_id_Ensembl <- ensembl_genes_from_seq_data[match(rownames(normalized_counts_dataframe), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_id_Ensembl
normalized_counts_dataframe$gene_id_Ensembl_remove_version_info <- ensembl_genes_from_seq_data[match(rownames(normalized_counts_dataframe), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$gene_id_Ensembl_remove_version_info
normalized_counts_dataframe$chr_name <- ensembl_genes_from_seq_data[match(rownames(normalized_counts_dataframe), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$chr_name
normalized_counts_dataframe$genomic_strand <- ensembl_genes_from_seq_data[match(rownames(normalized_counts_dataframe), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_strand
normalized_counts_dataframe$genomic_start <- ensembl_genes_from_seq_data[match(rownames(normalized_counts_dataframe), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_start
normalized_counts_dataframe$genomic_end <- ensembl_genes_from_seq_data[match(rownames(normalized_counts_dataframe), ensembl_genes_from_seq_data$Gene_name_in_seq_data), ]$genomic_end
write.csv(normalized_counts_dataframe, file = paste(baseDir, "/results/results_Gencode_DESeq2_normalized_counts.csv", sep=''))

###########################################################################################################################################

###########################################Inspect and check the data ####################################################################
# Take a look at the results of DESeq2. What do the columns contain?
mcols(res_all, use.name=T)
dim(res_all)


# Plot step 1: volcano plots ----------------------------------------------

#####################################################################################################################
                      # Plot Step 1: volcano plots
#####################################################################################################################
######diagnostic plot
# plotMA plot log2 fold change to basemean of normalized counts, and red dots mean genes with an adjusted p value below
# a threshold (here 0.1, the default). For other plot functions, check beginner_DESeq2.pres_all_dataframe_step1
plotMA(res_all, ylim = c(-10, 10))
# Volcano plots to visualize significant genes
res_all_dataframe_step1 <- data.frame(res_all)
range(res_all_dataframe_step1$log2FoldChange, na.rm = T)
range(-log10(res_all_dataframe_step1$padj), na.rm = T)
# since there are padj=0, which causes issues when plotting. So change those padj to a small number.
# see which genes have padj=0
res_all_dataframe_step1[which(res_all_dataframe_step1$padj==0), ] # 3 genes, Olfr446, Olfr536, and Olfr51 have padj = 0
# next check the minimum padj above 0
min(res_all_dataframe_step1[which(res_all_dataframe_step1$padj>0), ]$padj) # the minimum is 2.389498e-289
# next, set the 3 genes with padj=0 to 1e-300, otherwise it will cause issue when plotting with -log10(padj)
res_all_dataframe_step1[which(res_all_dataframe_step1$padj==0), ]$padj <- 1e-300
# confirm if they are correct
res_all_dataframe_step1[which(res_all_dataframe_step1$padj== 1e-300), ]

# check how many genes have threshold padj<0.05
sum(res_all$padj<0.05, na.rm = T)
sum(res_functional_OR$padj<0.05, na.rm = T)
sum(res_pseudogene_OR$padj<0.05, na.rm = T)
sum(res_TAAR$padj<0.05, na.rm = T)
sum(res_Ms4a$padj<0.05, na.rm = T)

plot_step1_1 <- ggplot(data=res_all_dataframe_step1, aes(x=log2FoldChange, y=-log10(padj), colour=res_all$padj<0.05)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_point(alpha=0.75, pch=16, size=1) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  geom_hline(yintercept=1.30103, linetype="dashed", 
             color = "black", size=0.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Trim66_KO homo vs het with threshold padj < 0.05, total: 2910; \nOR functional genes: 1025 of 1137, OR pseudogene genes: 38 of 271, \nTAAR functional genes: 11 of 14, no Ms4a genes') +
  theme(plot.title = element_text(vjust = 1, size = 5)) + # adjust the title size
  xlab("log2 fold change") + ylab("-log10 padj") + xlim(-8, 10) + ylim(-1,320) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
plot_step1_1

# add text to the data with signicant changes in the above figure
res_all_dataframe_step1$name <- rownames(res_all_dataframe_step1)
plot_step1_2 <- ggplot(data=res_all_dataframe_step1, aes(x=log2FoldChange, y=-log10(padj), label=name, colour=padj<0.05)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_point(alpha=0.75, pch=16, size=1) +
  geom_text(data=subset(res_all_dataframe_step1, padj<0.05), size =1, colour='blue') +  # label=name in the above aes() suggest the text to label
  geom_hline(yintercept=1.30103, linetype="dashed", 
             color = "black", size=0.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Trim66_KO homo vs het with threshold padj < 0.05, total: 2910; \nOR functional genes: 1025 of 1137, OR pseudogene genes: 38 of 271, \nTAAR functional genes: 11 of 14, no Ms4a genes') +
  theme(plot.title = element_text(vjust = 1, size = 5)) + # adjust the title size
  xlab("log2 fold change") + ylab("-log10 padj") + xlim(-8, 10) + ylim(-1,320) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
plot_step1_2

grid.arrange(plot_step1_1, plot_step1_2, ncol=2)  # put the two figures side by side

# ggpubr is a popular solution, particularly for those utlising ggplot2 for their plotting needs. 
# The main feature from this package is ggarrange() which will help arrange plots and allows further customisation in removing labels (for shared x or y axes), combining legends and even labelling each plot.
arrange_plots_step1_1and2 <- ggarrange(plot_step1_1, plot_step1_2, ncol = 2, nrow = 1)
ggsave(paste(baseDir, "/results/volcano_plots_padj_0.05.pdf", sep=''), arrange_plots_step1_1and2, width = 30, height = 15, units ="cm")

# check how many genes have threshold padj<0.05 and abs(FoldChange)>=1.5 
sum(res_all$padj<0.05 & abs(res_all$log2FoldChange)>=0.585, na.rm = T)
sum(res_functional_OR$padj<0.05 & abs(res_functional_OR$log2FoldChange)>=0.585, na.rm = T)
sum(res_pseudogene_OR$padj<0.05 & abs(res_pseudogene_OR$log2FoldChange)>=0.585, na.rm = T)
sum(res_TAAR$padj<0.05 & abs(res_TAAR$log2FoldChange)>=0.585, na.rm = T)
sum(res_Ms4a$padj<0.05 & abs(res_Ms4a$log2FoldChange)>=0.585, na.rm = T)

plot_step1_3 <- ggplot(data=res_all_dataframe_step1, aes(x=log2FoldChange, y=-log10(padj), colour=res_all$padj<0.05 & abs(res_all$log2FoldChange)>=0.585)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_point(alpha=0.75, pch=16, size=1) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  geom_hline(yintercept=1.30103, linetype="dashed", 
             color = "black", size=0.5) +
  geom_vline(xintercept=c(-0.585, 0.585), linetype="dashed",  # c(-0.585, 0.585) to add two vertical lines
             color = "black", size=0.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Trim66_KO homo vs het with threshold padj < 0.05 and abs(FoldChange)>= 1.5 \ntotal: 1384; \nOR functional genes: 1016 of 1137, OR pseudogene genes: 37 of 271, \nTAAR functional genes: 11 of 14, no Ms4a genes') +
  theme(plot.title = element_text(vjust = 1, size = 5)) + # adjust the title position and size
  xlab("log2 fold change") + ylab("-log10 padj") + xlim(-8, 10) + ylim(-1,320) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
plot_step1_3

# add text to the data with signicant changes in the above figure
res_all_dataframe_step1$name <- rownames(res_all_dataframe_step1)
plot_step1_4 <- ggplot(data=res_all_dataframe_step1, aes(x=log2FoldChange, y=-log10(padj), label=name, colour=padj<0.05 & abs(log2FoldChange)>=0.585)) +
  scale_color_manual(values = c("grey", "red")) +
  geom_point(alpha=0.75, pch=16, size=1) +
  geom_text(data=subset(res_all_dataframe_step1, (padj<0.05 & abs(log2FoldChange)>=0.585) ==TRUE), size =1, colour='blue') +  # label=name in the above aes() suggest the text to label
  geom_hline(yintercept=1.30103, linetype="dashed", 
             color = "black", size=0.5) +
  geom_vline(xintercept=c(-0.585, 0.585), linetype="dashed",  # c(-0.585, 0.585) to add two vertical lines
             color = "black", size=0.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Trim66_KO homo vs het with threshold padj < 0.05 and abs(FoldChange)>= 1.5 \ntotal: 1384; \nOR functional genes: 1016 of 1137, OR pseudogene genes: 37 of 271, \nTAAR functional genes: 11 of 14, no Ms4a genes') +
  theme(plot.title = element_text(vjust = 1, size = 5)) + # adjust the title position and size
  xlab("log2 fold change") + ylab("-log10 padj") + xlim(-8, 10) + ylim(-1,320) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
plot_step1_4

grid.arrange(plot_step1_3, plot_step1_4, ncol=2)  # put the two figures side by side

# ggpubr is a popular solution, particularly for those utlising ggplot2 for their plotting needs. 
# The main feature from this package is ggarrange() which will help arrange plots and allows further customisation in removing labels (for shared x or y axes), combining legends and even labelling each plot.
arrange_plots_step1_3and4 <- ggarrange(plot_step1_3, plot_step1_4, ncol = 2, nrow = 1)
ggsave(paste(baseDir, "/results/volcano_plots_padj_0.05_foldChange_1.5.pdf", sep=''), arrange_plots_step1_3and4, width = 30, height = 15, units ="cm")

# highlight functional ORs (red1), pseudogene ORs (springgreen3), functional TAARs (blue2),and pseudogene TAARs (magenta1)
res_all_dataframe_step1_5and6 <- res_all_dataframe_step1 %>% 
  mutate(is_OR_TAAR_pseudogene = "non_olfactory_receptor") %>% 
  mutate(is_OR_TAAR_pseudogene = replace(is_OR_TAAR_pseudogene, grepl("protein_coding", gene_type) & grepl("Olfr", name), "functional_OR")) %>% 
  mutate(is_OR_TAAR_pseudogene = replace(is_OR_TAAR_pseudogene, grepl("pseudogene", gene_type) & grepl("Olfr", name), "pseudogene_OR")) %>% 
  mutate(is_OR_TAAR_pseudogene = replace(is_OR_TAAR_pseudogene, grepl("protein_coding", gene_type) & grepl("Taar[2-9]", name), "functional_TAAR")) %>% 
  mutate(is_OR_TAAR_pseudogene = replace(is_OR_TAAR_pseudogene, grepl("pseudogene", gene_type) & grepl("Taar[2-9]", name), "pseudogene_TAAR"))
plot_step1_5 <- ggplot(data=res_all_dataframe_step1_5and6, aes(x=log2FoldChange, y=-log10(padj), colour=is_OR_TAAR_pseudogene, alpha=is_OR_TAAR_pseudogene)) +
  geom_point(pch=16, size=0.5) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_color_manual(values = c("red1", "blue2", "grey40", "springgreen3", "magenta1")) +
  scale_alpha_manual(values = c(1, 1, 0.5, 1, 1)) +
  geom_hline(yintercept=1.30103, linetype="dashed", 
             color = "black", size=0.5) +
  geom_vline(xintercept=c(-0.585, 0.585), linetype="dashed",  # c(-0.585, 0.585) to add two vertical lines
             color = "black", size=0.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Trim66_KO homo vs het with receptors highlighted: functional ORs (red1), pseudogene ORs (springgreen3), functional TAARs (blue2),and pseudogene TAARs (magenta1)\nthreshold padj < 0.05 and abs(FoldChange)>= 1.5 \ntotal: 1384; \nOR functional genes: 1016 of 1137, OR pseudogene genes: 37 of 271, \nTAAR functional genes: 11 of 14, no Ms4a genes') +
  theme(plot.title = element_text(vjust = 1, size = 5)) + # adjust the title position and size
  xlab("log2 fold change") + ylab("-log10 padj") + xlim(-8, 10) + ylim(-1,320) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
plot_step1_5

# add text to the data with receptor names in the above figure
plot_step1_6 <- ggplot(data=res_all_dataframe_step1_5and6, aes(x=log2FoldChange, y=-log10(padj), colour=is_OR_TAAR_pseudogene, alpha=is_OR_TAAR_pseudogene, label=name)) +
  geom_point(pch=16, size=0.5) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_color_manual(values = c("red1", "blue2", "grey40", "springgreen3", "magenta1")) +
  scale_alpha_manual(values = c(1, 1, 0.5, 1, 1)) +
  geom_text(data=subset(res_all_dataframe_step1_5and6, !grepl("non_olfactory_receptor", is_OR_TAAR_pseudogene)), size =1, colour='black') +  # label=name in the above aes() suggest the text to label
  geom_hline(yintercept=1.30103, linetype="dashed", 
             color = "black", size=0.5) +
  geom_vline(xintercept=c(-0.585, 0.585), linetype="dashed",  # c(-0.585, 0.585) to add two vertical lines
             color = "black", size=0.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(0.75)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Trim66_KO homo vs het with receptors highlighted: functional ORs (red1), pseudogene ORs (springgreen3), functional TAARs (blue2),and pseudogene TAARs (magenta1)\nthreshold padj < 0.05 and abs(FoldChange)>= 1.5 \ntotal: 1384; \nOR functional genes: 1016 of 1137, OR pseudogene genes: 37 of 271, \nTAAR functional genes: 11 of 14, no Ms4a genes') +
  theme(plot.title = element_text(vjust = 1, size = 5)) + # adjust the title position and size
  xlab("log2 fold change") + ylab("-log10 padj") + xlim(-8, 10) + ylim(-1,320) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
plot_step1_6

grid.arrange(plot_step1_5, plot_step1_6, ncol=2)  # put the two figures side by side

# ggpubr is a popular solution, particularly for those utlising ggplot2 for their plotting needs. 
# The main feature from this package is ggarrange() which will help arrange plots and allows further customisation in removing labels (for shared x or y axes), combining legends and even labelling each plot.
arrange_plots_step1_5and6 <- ggarrange(plot_step1_5, plot_step1_6, ncol = 2, nrow = 1)
ggsave(paste(baseDir, "/results/volcano_plots_olfactory_receptors_highlighted.pdf", sep=''), arrange_plots_step1_5and6, width = 30, height = 15, units ="cm")

# Plot step 2: pie charts and stacked bar plots ---------------------------

#####################################################################################################################
                          # Plot Step 2: pie chart and stacked bar plot for percentage of downregulated and upregulated receptor genes
     # include all functional OR genes, class I OR, class II OR, TAAR, and OR pseudogenes
#####################################################################################################################
all_genes_number_step2 <- nrow(res_all)
# upregualted genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_all_genes_number_step2 <- sum(res_all$padj<0.05 & res_all$log2FoldChange>=0.585, na.rm = T)
# downregualted genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_all_genes_number_step2 <- sum(res_all$padj<0.05 & res_all$log2FoldChange<=-0.585, na.rm = T)
# other genes are counted as not changed
nochange_all_genes_number_step2 <- all_genes_number_step2-upregulated_all_genes_number_step2-downregulated_all_genes_number_step2

# functional OR genes
all_functional_OR_genes_number_step2 <- nrow(res_functional_OR)
# upregualted functional OR genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_functional_OR_genes_number_step2 <- sum(res_functional_OR$padj<0.05 & res_functional_OR$log2FoldChange>=0.585, na.rm = T)
# downregualted functional OR genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_functional_OR_genes_number_step2 <- sum(res_functional_OR$padj<0.05 & res_functional_OR$log2FoldChange<=-0.585, na.rm = T)
# other functional OR genes are counted as not changed
nochange_functional_OR_genes_number_step2 <- all_functional_OR_genes_number_step2-upregulated_functional_OR_genes_number_step2-downregulated_functional_OR_genes_number_step2

# pseudogene OR genes
all_pseudogene_OR_genes_number_step2 <- nrow(res_pseudogene_OR)
# upregualted pseudogene OR genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_pseudogene_OR_genes_number_step2 <- sum(res_pseudogene_OR$padj<0.05 & res_pseudogene_OR$log2FoldChange>=0.585, na.rm = T)
# downregualted pseudogene OR genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_pseudogene_OR_genes_number_step2 <- sum(res_pseudogene_OR$padj<0.05 & res_pseudogene_OR$log2FoldChange<=-0.585, na.rm = T)
# other pseudogene OR genes are counted as not changed
nochange_pseudogene_OR_genes_number_step2 <- all_pseudogene_OR_genes_number_step2-upregulated_pseudogene_OR_genes_number_step2-downregulated_pseudogene_OR_genes_number_step2

# functional TAAR genes. Omit TAAR1 and TAAR7c.ps. These 2 genes are not significantly changed here.
all_functional_TAAR_genes_number_step2 <- 14
# upregualted functional TAAR genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_functional_TAAR_genes_number_step2 <- sum(!grepl("^Taar[17]\\.", rownames(res_TAAR)) & 
                                                        res_TAAR$padj<0.05 & 
                                                        res_TAAR$log2FoldChange>=0.585, na.rm = T)
# downregualted functional TAAR genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_functional_TAAR_genes_number_step2 <- sum(!grepl("^Taar[17]\\.", rownames(res_TAAR)) &
                                                          res_TAAR$padj<0.05 & 
                                                          res_TAAR$log2FoldChange<=-0.585, na.rm = T)
# other functional TAAR genes are counted as not changed
nochange_functional_TAAR_genes_number_step2 <- all_functional_TAAR_genes_number_step2-upregulated_functional_TAAR_genes_number_step2-downregulated_functional_TAAR_genes_number_step2

#####################################################################################################################
            # check class I and class II OR genes
#####################################################################################################################
# load class I OR genes from Longzhi's Chemical Senses paper in 2018, note that he used GRCm38.p5
classI_OR_list <- read.delim("/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/ClassI_ORs_Longzhi_paper_2018.csv", header=T, sep=",",
                             row.names=NULL, as.is=T)
all(classI_OR_list$gene.id.version.in.GRCm38.p6 %in% res_OR$gene_id_Ensembl)
# note that there is one gene olfr1527-ps1 is not named OR in GRCm38.p6. I am not sure if it is OR, as in Junji's NC paper, they claimed 158 class I genes, but Longzhi claimed 159 genes. Maybe this is the extra one.
all(classI_OR_list$gene.name %in% rownames(res_OR))
# remove olfr1527-ps1 from class I OR list
classI_OR_list_remove_Olfr1527.ps1 <- classI_OR_list[!grepl("Olfr1527", classI_OR_list$gene.name), ]

# how many class I OR genes show signicant changes. Total 158 class I OR genes.
res_ClassI_OR <- res_OR[res_OR$gene_id_Ensembl %in% classI_OR_list$gene.id.version.in.GRCm38.p6, ]
dim(res_ClassI_OR)
dim(classI_OR_list_remove_Olfr1527.ps1)
write.csv(res_ClassI_OR, file = paste(baseDir, "/results/results_ClassI_OR_genes_Gencode_Trim66_homo_vs_het.csv", sep=''))

# functional ClassI_OR genes
all_functional_ClassI_OR_genes_number_step2 <- nrow(res_ClassI_OR[grepl("protein_coding", res_ClassI_OR$gene_type), ])
# upregualted functional ClassI_OR genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_functional_ClassI_OR_genes_number_step2 <- sum(res_ClassI_OR$padj<0.05 & res_ClassI_OR$log2FoldChange>=0.585 & grepl("protein_coding", res_ClassI_OR$gene_type), na.rm = T)
# downregualted functional ClassI_OR genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_functional_ClassI_OR_genes_number_step2 <- sum(res_ClassI_OR$padj<0.05 & res_ClassI_OR$log2FoldChange<=-0.585 & grepl("protein_coding", res_ClassI_OR$gene_type), na.rm = T)
# other functional ClassI_OR genes are counted as not changed
nochange_functional_ClassI_OR_genes_number_step2 <- all_functional_ClassI_OR_genes_number_step2-upregulated_functional_ClassI_OR_genes_number_step2-downregulated_functional_ClassI_OR_genes_number_step2

# pseudogene ClassI_OR genes
all_pseudogene_ClassI_OR_genes_number_step2 <- nrow(res_ClassI_OR[grepl("pseudogene", res_ClassI_OR$gene_type), ])
# upregualted pseudogene ClassI_OR genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_pseudogene_ClassI_OR_genes_number_step2 <- sum(res_ClassI_OR$padj<0.05 & res_ClassI_OR$log2FoldChange>=0.585 & grepl("pseudogene", res_ClassI_OR$gene_type), na.rm = T)
# downregualted pseudogene ClassI_OR genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_pseudogene_ClassI_OR_genes_number_step2 <- sum(res_ClassI_OR$padj<0.05 & res_ClassI_OR$log2FoldChange<=-0.585 & grepl("pseudogene", res_ClassI_OR$gene_type), na.rm = T)
# other pseudogene ClassI_OR genes are counted as not changed
nochange_pseudogene_ClassI_OR_genes_number_step2 <- all_pseudogene_ClassI_OR_genes_number_step2-upregulated_pseudogene_ClassI_OR_genes_number_step2-downregulated_pseudogene_ClassI_OR_genes_number_step2

#######################################################################
# I tried to load Longzhi's data that have 1243 matches with res_OR. 
# As it is supposed as 1408-158=1250 class II genes, so I will directly use the remaining OR genes as class II genes
#######################################################################
res_ClassII_OR <- res_OR[!res_OR$gene_id_Ensembl %in% classI_OR_list$gene.id.version.in.GRCm38.p6, ]
dim(res_ClassII_OR)
write.csv(res_ClassII_OR, file = paste(baseDir, "/results/results_ClassII_OR_genes_Gencode_Trim66_homo_vs_het.csv", sep=''))

# functional ClassII_OR genes
all_functional_ClassII_OR_genes_number_step2 <- nrow(res_ClassII_OR[grepl("protein_coding", res_ClassII_OR$gene_type), ])
# upregualted functional ClassII_OR genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_functional_ClassII_OR_genes_number_step2 <- sum(res_ClassII_OR$padj<0.05 & res_ClassII_OR$log2FoldChange>=0.585 & grepl("protein_coding", res_ClassII_OR$gene_type), na.rm = T)
# downregualted functional ClassII_OR genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_functional_ClassII_OR_genes_number_step2 <- sum(res_ClassII_OR$padj<0.05 & res_ClassII_OR$log2FoldChange<=-0.585 & grepl("protein_coding", res_ClassII_OR$gene_type), na.rm = T)
# other functional ClassII_OR genes are counted as not changed
nochange_functional_ClassII_OR_genes_number_step2 <- all_functional_ClassII_OR_genes_number_step2-upregulated_functional_ClassII_OR_genes_number_step2-downregulated_functional_ClassII_OR_genes_number_step2

# pseudogene ClassII_OR genes
all_pseudogene_ClassII_OR_genes_number_step2 <- nrow(res_ClassII_OR[grepl("pseudogene", res_ClassII_OR$gene_type), ])
# upregualted pseudogene ClassII_OR genes with criteria of padj<0.05 and foldchange>=1.5
upregulated_pseudogene_ClassII_OR_genes_number_step2 <- sum(res_ClassII_OR$padj<0.05 & res_ClassII_OR$log2FoldChange>=0.585 & grepl("pseudogene", res_ClassII_OR$gene_type), na.rm = T)
# downregualted pseudogene ClassII_OR genes with criteria of padj<0.05 and foldchange<=0.67
downregulated_pseudogene_ClassII_OR_genes_number_step2 <- sum(res_ClassII_OR$padj<0.05 & res_ClassII_OR$log2FoldChange<=-0.585 & grepl("pseudogene", res_ClassII_OR$gene_type), na.rm = T)
# other pseudogene ClassII_OR genes are counted as not changed
nochange_pseudogene_ClassII_OR_genes_number_step2 <- all_pseudogene_ClassII_OR_genes_number_step2-upregulated_pseudogene_ClassII_OR_genes_number_step2-downregulated_pseudogene_ClassII_OR_genes_number_step2


#######################################################################
Note: For record.  Do not run!
#######################################################################
#  Do not run this. This is just for record.
# I tried to load Longzhi's data that have 1243 matches with res_OR. 
# As it is supposed as 1408-158=1250 class II genes, so I will directly use the remaining OR genes as class II genes
#      Skip this part!!!!!!!!!!!!!
#######################################################################
# load class II OR genes from Longzhi's Chemical Senses paper in 2018
classII_OR_list <- read.delim("/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/ClassII_ORs_Longzhi_paper_2018.csv", header=T, sep=",",
                              row.names=NULL, as.is=T)
all(classII_OR_list$gene.id.version.in.GRCm38.p6 %in% res_OR$gene_id_Ensembl)
# note that there are 4 genes including Olfr1010, Olfr892-ps1, Olfr964-ps1, and Olfr192 that have two version in GRCm38.p5 but only one in GRCm38.p6
# and Olfr18 and Olfr912 have two versions in GRCm38.p5 but correspond to differnt names in GRCm38.p6. I generated this new file by matching the gene ID in GRCm38.p5 to GRCm38.p6, then I added the new gene ID version and new gene name in GRCm38.p6 in Longzhi's original list.
classII_OR_list_corrected <- classII_OR_list[classII_OR_list$gene.id.version.in.GRCm38.p6 %in% res_OR$gene_id_Ensembl, ]
dim(classII_OR_list_corrected)  #Total 1243 class II OR genes.
#######################################################################
              End of the note!
#######################################################################


#####################################################################################################################
                         # check if zones of OR genes matter
#####################################################################################################################
# load class II OR file from Longzhi's Chemical Senses paper in 2018 containing the zone information
classII_OR_list <- read.delim("/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/ClassII_ORs_Longzhi_paper_2018.csv", header=T, sep=",",
                              row.names=NULL, as.is=T)
all(classII_OR_list$gene.id.version.in.GRCm38.p6 %in% res_OR$gene_id_Ensembl)
#################################################################################################
# note that there are 4 genes including Olfr1010, Olfr892-ps1, Olfr964-ps1, and Olfr192 that have two version in GRCm38.p5 but only one in GRCm38.p6
# and Olfr18 and Olfr912 have two versions in GRCm38.p5 but correspond to differnt names in GRCm38.p6. I generated this new file by matching the gene ID in GRCm38.p5 to GRCm38.p6, then I added the new gene ID version and new gene name in GRCm38.p6 in Longzhi's original list.
classII_OR_list_corrected <- classII_OR_list[classII_OR_list$gene.id.version.in.GRCm38.p6 %in% res_OR$gene_id_Ensembl, ]
dim(classII_OR_list_corrected)  #Total 1243 class II OR genes.
#################################################################################################
# add one column in res_all to indicate Class I/II
res_all$OR_Class <- NA
res_all[res_all$gene_id_Ensembl %in% res_ClassI_OR$gene_id_Ensembl, ]$OR_Class <- "Class_I"
res_all[which(res_all$OR_Class == "Class_I"), ]
sum(res_all$OR_Class == "Class_I", na.rm = T)
res_all[res_all$gene_id_Ensembl %in% res_ClassII_OR$gene_id_Ensembl, ]$OR_Class <- "Class_II"
sum(res_all$OR_Class == "Class_II", na.rm = T)
# add one column in res_all to indicate zone information
res_all$inferred_zone_index_Longzhi_2018 <- NA
all(classI_OR_list$gene.id.version.in.GRCm38.p6 %in% res_all$gene_id_Ensembl) #all class I ORs can be found in res_all
res_all[match(classI_OR_list$gene.id.version.in.GRCm38.p6, res_all$gene_id_Ensembl), ]$inferred_zone_index_Longzhi_2018 <- classI_OR_list$inferred.zone.index
all(classII_OR_list$gene.id.version.in.GRCm38.p6 %in% res_all$gene_id_Ensembl) #not all class II ORs can be found in res_all
sum(classII_OR_list$gene.id.version.in.GRCm38.p6 %in% res_all$gene_id_Ensembl) # 1253 matches
dim(classII_OR_list) # there are 1258 class II ORs in Longzhi's list
# remove the unmatched ones in the lis
classII_OR_list_step2 <- classII_OR_list
classII_OR_list_step2 <- classII_OR_list_step2[classII_OR_list_step2$gene.id.version.in.GRCm38.p6 %in% res_all$gene_id_Ensembl, ]
all(classII_OR_list_step2$gene.id.version.in.GRCm38.p6 %in% res_all$gene_id_Ensembl) # now all matched
sum(classII_OR_list_step2$gene.id.version.in.GRCm38.p6 %in% res_all$gene_id_Ensembl)
res_all[match(classII_OR_list_step2$gene.id.version.in.GRCm38.p6, res_all$gene_id_Ensembl), ]$inferred_zone_index_Longzhi_2018 <- classII_OR_list_step2$inferred.zone.index
# add one column to use zone 1-5 instead of zone index, note that 4.5 will be round to 4 instead of 5 in R
res_all$zone_1to5_by_Qian_rounded_in_R <- round(as.numeric(res_all$inferred_zone_index_Longzhi_2018))
sum(res_all$zone_1to5_by_Qian_rounded_in_R >=1, na.rm = T)
# save to file
write.csv(res_all, 
          file = paste(baseDir, 
                       "/results/results_Gencode_Trim66_homo_vs_het_with_Class_and_zone_info.csv", 
                       sep=''))
# how many ORs have the zone index
res_all_ORs_with_zone_info_step2 <- res_all[grepl("^[1-5]",res_all$inferred_zone_index_Longzhi_2018), ]
dim(res_all_ORs_with_zone_info_step2) #1011 ORs have zone information
sum(res_all_ORs_with_zone_info_step2$OR_Class == "Class_I", na.rm = T)  #115 Class I ORs have zone information

# how many ORs in different zones are significantly changed
unique(res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R)
for (i in 1:5) {
  assign(paste("all_functional_zone_", i, "_OR_genes_number_step2", sep = ""), value = nrow(res_all_ORs_with_zone_info_step2[grepl("protein_coding", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, ]))
  assign(paste("upregulated_functional_zone_", i, "_OR_genes_number_step2", sep = ""), value = sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange>=0.585 & grepl("protein_coding", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T))
  assign(paste("downregulated_functional_zone_", i, "_OR_genes_number_step2", sep = ""), value = sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange<=-0.585 & grepl("protein_coding", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T))
  assign(paste("nochange_functional_zone_", i, "_OR_genes_number_step2", sep = ""), 
         value=nrow(res_all_ORs_with_zone_info_step2[grepl("protein_coding", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, ]) -
           sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange>=0.585 & grepl("protein_coding", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T) -
           sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange<=-0.585 & grepl("protein_coding", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T))
  assign(paste("all_pseudogene_zone_", i, "_OR_genes_number_step2", sep = ""), value = nrow(res_all_ORs_with_zone_info_step2[grepl("pseudogene", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, ]))
  assign(paste("upregulated_pseudogene_zone_", i, "_OR_genes_number_step2", sep = ""), value = sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange>=0.585 & grepl("pseudogene", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T))
  assign(paste("downregulated_pseudogene_zone_", i, "_OR_genes_number_step2", sep = ""), value = sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange<=-0.585 & grepl("pseudogene", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T))
  assign(paste("nochange_pseudogene_zone_", i, "_OR_genes_number_step2", sep = ""), 
         value = nrow(res_all_ORs_with_zone_info_step2[grepl("pseudogene", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, ]) -
           sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange>=0.585 & grepl("pseudogene", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T) - 
           sum(res_all_ORs_with_zone_info_step2$padj<0.05 & res_all_ORs_with_zone_info_step2$log2FoldChange<=-0.585 & grepl("pseudogene", res_all_ORs_with_zone_info_step2$gene_type) & res_all_ORs_with_zone_info_step2$zone_1to5_by_Qian_rounded_in_R == i, na.rm = T))
}

################ Start  Pie chart ######
# Create data for the graph.
all_genes_pie_data_step2 <-  c(upregulated_all_genes_number_step2, downregulated_all_genes_number_step2, nochange_all_genes_number_step2)
pie_labels_step2 <-  c("Upregulated genes (padj<0.05 and FoldChange>=1.5)","Downregulated genes (padj<0.05 and FoldChange<=-1.5)","Other genes")
all_genes_pie_data_percent_step2<- paste(round(100*all_genes_pie_data_step2/sum(all_genes_pie_data_step2), 1), 
                                         "%, ", 
                                         all_genes_pie_data_step2, 
                                         " genes",
                                         sep = "")

all_functional_OR_genes_pie_data_step2 <-  c(upregulated_functional_OR_genes_number_step2, downregulated_functional_OR_genes_number_step2, nochange_functional_OR_genes_number_step2)
all_functional_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_OR_genes_pie_data_step2/sum(all_functional_OR_genes_pie_data_step2), 1), 
                                                       "%, ", 
                                                       all_functional_OR_genes_pie_data_step2, 
                                                       " genes",
                                                       sep = "")

all_pseudogene_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_OR_genes_number_step2, downregulated_pseudogene_OR_genes_number_step2, nochange_pseudogene_OR_genes_number_step2)
all_pseudogene_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_OR_genes_pie_data_step2/sum(all_pseudogene_OR_genes_pie_data_step2), 1), 
                                                       "%, ", 
                                                       all_pseudogene_OR_genes_pie_data_step2, 
                                                       " genes",
                                                       sep = "")

all_functional_ClassI_OR_genes_pie_data_step2 <-  c(upregulated_functional_ClassI_OR_genes_number_step2, downregulated_functional_ClassI_OR_genes_number_step2, nochange_functional_ClassI_OR_genes_number_step2)
all_functional_ClassI_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_ClassI_OR_genes_pie_data_step2/sum(all_functional_ClassI_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_functional_ClassI_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_pseudogene_ClassI_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_ClassI_OR_genes_number_step2, downregulated_pseudogene_ClassI_OR_genes_number_step2, nochange_pseudogene_ClassI_OR_genes_number_step2)
all_pseudogene_ClassI_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_ClassI_OR_genes_pie_data_step2/sum(all_pseudogene_ClassI_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_pseudogene_ClassI_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_functional_ClassII_OR_genes_pie_data_step2 <-  c(upregulated_functional_ClassII_OR_genes_number_step2, downregulated_functional_ClassII_OR_genes_number_step2, nochange_functional_ClassII_OR_genes_number_step2)
all_functional_ClassII_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_ClassII_OR_genes_pie_data_step2/sum(all_functional_ClassII_OR_genes_pie_data_step2), 1), 
                                                               "%, ", 
                                                               all_functional_ClassII_OR_genes_pie_data_step2, 
                                                               " genes",
                                                               sep = "")

all_pseudogene_ClassII_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_ClassII_OR_genes_number_step2, downregulated_pseudogene_ClassII_OR_genes_number_step2, nochange_pseudogene_ClassII_OR_genes_number_step2)
all_pseudogene_ClassII_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_ClassII_OR_genes_pie_data_step2/sum(all_pseudogene_ClassII_OR_genes_pie_data_step2), 1), 
                                                               "%, ", 
                                                               all_pseudogene_ClassII_OR_genes_pie_data_step2, 
                                                               " genes",
                                                               sep = "")

all_functional_TAAR_genes_pie_data_step2 <-  c(upregulated_functional_TAAR_genes_number_step2, downregulated_functional_TAAR_genes_number_step2, nochange_functional_TAAR_genes_number_step2)
all_functional_TAAR_genes_pie_data_percent_step2<- paste(round(100*all_functional_TAAR_genes_pie_data_step2/sum(all_functional_TAAR_genes_pie_data_step2), 1), 
                                                         "%, ", 
                                                         all_functional_TAAR_genes_pie_data_step2, 
                                                         " genes",
                                                         sep = "")

all_functional_zone_1_OR_genes_pie_data_step2 <-  c(upregulated_functional_zone_1_OR_genes_number_step2, downregulated_functional_zone_1_OR_genes_number_step2, nochange_functional_zone_1_OR_genes_number_step2)
all_functional_zone_1_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_zone_1_OR_genes_pie_data_step2/sum(all_functional_zone_1_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_functional_zone_1_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_pseudogene_zone_1_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_zone_1_OR_genes_number_step2, downregulated_pseudogene_zone_1_OR_genes_number_step2, nochange_pseudogene_zone_1_OR_genes_number_step2)
all_pseudogene_zone_1_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_zone_1_OR_genes_pie_data_step2/sum(all_pseudogene_zone_1_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_pseudogene_zone_1_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_functional_zone_2_OR_genes_pie_data_step2 <-  c(upregulated_functional_zone_2_OR_genes_number_step2, downregulated_functional_zone_2_OR_genes_number_step2, nochange_functional_zone_2_OR_genes_number_step2)
all_functional_zone_2_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_zone_2_OR_genes_pie_data_step2/sum(all_functional_zone_2_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_functional_zone_2_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_pseudogene_zone_2_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_zone_2_OR_genes_number_step2, downregulated_pseudogene_zone_2_OR_genes_number_step2, nochange_pseudogene_zone_2_OR_genes_number_step2)
all_pseudogene_zone_2_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_zone_2_OR_genes_pie_data_step2/sum(all_pseudogene_zone_2_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_pseudogene_zone_2_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_functional_zone_3_OR_genes_pie_data_step2 <-  c(upregulated_functional_zone_3_OR_genes_number_step2, downregulated_functional_zone_3_OR_genes_number_step2, nochange_functional_zone_3_OR_genes_number_step2)
all_functional_zone_3_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_zone_3_OR_genes_pie_data_step2/sum(all_functional_zone_3_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_functional_zone_3_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_pseudogene_zone_3_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_zone_3_OR_genes_number_step2, downregulated_pseudogene_zone_3_OR_genes_number_step2, nochange_pseudogene_zone_3_OR_genes_number_step2)
all_pseudogene_zone_3_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_zone_3_OR_genes_pie_data_step2/sum(all_pseudogene_zone_3_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_pseudogene_zone_3_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_functional_zone_4_OR_genes_pie_data_step2 <-  c(upregulated_functional_zone_4_OR_genes_number_step2, downregulated_functional_zone_4_OR_genes_number_step2, nochange_functional_zone_4_OR_genes_number_step2)
all_functional_zone_4_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_zone_4_OR_genes_pie_data_step2/sum(all_functional_zone_4_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_functional_zone_4_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_pseudogene_zone_4_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_zone_4_OR_genes_number_step2, downregulated_pseudogene_zone_4_OR_genes_number_step2, nochange_pseudogene_zone_4_OR_genes_number_step2)
all_pseudogene_zone_4_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_zone_4_OR_genes_pie_data_step2/sum(all_pseudogene_zone_4_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_pseudogene_zone_4_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_functional_zone_5_OR_genes_pie_data_step2 <-  c(upregulated_functional_zone_5_OR_genes_number_step2, downregulated_functional_zone_5_OR_genes_number_step2, nochange_functional_zone_5_OR_genes_number_step2)
all_functional_zone_5_OR_genes_pie_data_percent_step2<- paste(round(100*all_functional_zone_5_OR_genes_pie_data_step2/sum(all_functional_zone_5_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_functional_zone_5_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

all_pseudogene_zone_5_OR_genes_pie_data_step2 <-  c(upregulated_pseudogene_zone_5_OR_genes_number_step2, downregulated_pseudogene_zone_5_OR_genes_number_step2, nochange_pseudogene_zone_5_OR_genes_number_step2)
all_pseudogene_zone_5_OR_genes_pie_data_percent_step2<- paste(round(100*all_pseudogene_zone_5_OR_genes_pie_data_step2/sum(all_pseudogene_zone_5_OR_genes_pie_data_step2), 1), 
                                                              "%, ", 
                                                              all_pseudogene_zone_5_OR_genes_pie_data_step2, 
                                                              " genes",
                                                              sep = "")

# Give the chart file a name.
pdf(file = paste(resultDir,"/changed_genes_percentage_Trim66_KO_homo_vs_het.pdf", sep=''), paper = 'a4')
opar <- par() #save the original par setting
# note that if par(mfrow=c(4,2)) before pdf() will print 8 pie charts in different pages of pdf. 
# put par(mfrow=c(4,2)) after pdf() here will print 8 pie charts in the same page.
# Check https://www.statology.org/r-save-multiple-plots-to-pdf/
par(mfrow=c(4,2), mar=c(0,0,1,0)) # there are 8 pie charts, arrange them in one page with 3X3 by row. mar to set up the margins of figures.
# Plot the chart.
pie_col <- c("lightcoral", "springgreen3", "gray80")
pie(all_genes_pie_data_step2, 
    labels = all_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het all ", all_genes_number_step2, " genes", sep = ""),
    col = pie_col,
    cex.main=1)
legend("topleft", pie_labels_step2, cex = 0.5, fill = pie_col)

pie(all_functional_OR_genes_pie_data_step2, 
    labels = all_functional_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_OR_genes_number_step2, " functional OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_pseudogene_OR_genes_pie_data_step2, 
    labels = all_pseudogene_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_OR_genes_number_step2, " pseudogene OR genes", sep=""),
    col = pie_col,
    cex.main=1)

pie(all_functional_ClassI_OR_genes_pie_data_step2, 
    labels = all_functional_ClassI_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_ClassI_OR_genes_number_step2, " functional ClassI_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_pseudogene_ClassI_OR_genes_pie_data_step2, 
    labels = all_pseudogene_ClassI_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_ClassI_OR_genes_number_step2, " pseudogene ClassI_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_functional_ClassII_OR_genes_pie_data_step2, 
    labels = all_functional_ClassII_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_ClassII_OR_genes_number_step2, " functional ClassII_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_pseudogene_ClassII_OR_genes_pie_data_step2, 
    labels = all_pseudogene_ClassII_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_ClassII_OR_genes_number_step2, " pseudogene ClassII_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_functional_TAAR_genes_pie_data_step2, 
    labels = all_functional_TAAR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_TAAR_genes_number_step2, " functional TAAR genes", sep = ""),
    col = pie_col,
    cex.main=1)

#Plot another page of 10 plots of zone ORs
par(mfrow=c(5,2), mar=c(0,0,1,0))
pie(all_functional_zone_1_OR_genes_pie_data_step2, 
    labels = all_functional_zone_1_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_zone_1_OR_genes_number_step2, " functional zone_1_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)
legend("topleft", pie_labels_step2, cex = 0.5, fill = pie_col)

pie(all_pseudogene_zone_1_OR_genes_pie_data_step2, 
    labels = all_pseudogene_zone_1_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_zone_1_OR_genes_number_step2, " pseudogene zone_1_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_functional_zone_2_OR_genes_pie_data_step2, 
    labels = all_functional_zone_2_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_zone_2_OR_genes_number_step2, " functional zone_2_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_pseudogene_zone_2_OR_genes_pie_data_step2, 
    labels = all_pseudogene_zone_2_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_zone_2_OR_genes_number_step2, " pseudogene zone_2_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_functional_zone_3_OR_genes_pie_data_step2, 
    labels = all_functional_zone_3_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_zone_3_OR_genes_number_step2, " functional zone_3_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_pseudogene_zone_3_OR_genes_pie_data_step2, 
    labels = all_pseudogene_zone_3_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_zone_3_OR_genes_number_step2, " pseudogene zone_3_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_functional_zone_4_OR_genes_pie_data_step2, 
    labels = all_functional_zone_4_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_zone_4_OR_genes_number_step2, " functional zone_4_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_pseudogene_zone_4_OR_genes_pie_data_step2, 
    labels = all_pseudogene_zone_4_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_zone_4_OR_genes_number_step2, " pseudogene zone_4_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_functional_zone_5_OR_genes_pie_data_step2, 
    labels = all_functional_zone_5_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_functional_zone_5_OR_genes_number_step2, " functional zone_5_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

pie(all_pseudogene_zone_5_OR_genes_pie_data_step2, 
    labels = all_pseudogene_zone_5_OR_genes_pie_data_percent_step2, 
    main = paste("Trim66_KO_homo_vs_het ", all_pseudogene_zone_5_OR_genes_number_step2, " pseudogene zone_5_OR genes", sep = ""),
    col = pie_col,
    cex.main=1)

# Save the file.
dev.off()
par(opar) # back to the original par for next plotting

################ stacked bar plot by ggplot2 ######
all_data_bar_plot_step2 <- data.frame(Number_of_genes=c(all_genes_pie_data_step2, 
                                                        all_functional_TAAR_genes_pie_data_step2, 
                                                        all_functional_OR_genes_pie_data_step2, 
                                                        all_functional_ClassI_OR_genes_pie_data_step2, 
                                                        all_functional_ClassII_OR_genes_pie_data_step2, 
                                                        all_functional_zone_1_OR_genes_pie_data_step2, 
                                                        all_functional_zone_2_OR_genes_pie_data_step2, 
                                                        all_functional_zone_3_OR_genes_pie_data_step2, 
                                                        all_functional_zone_4_OR_genes_pie_data_step2, 
                                                        all_functional_zone_5_OR_genes_pie_data_step2, 
                                                        all_pseudogene_OR_genes_pie_data_step2, 
                                                        all_pseudogene_ClassI_OR_genes_pie_data_step2,
                                                        all_pseudogene_ClassII_OR_genes_pie_data_step2,
                                                        all_pseudogene_zone_1_OR_genes_pie_data_step2, 
                                                        all_pseudogene_zone_2_OR_genes_pie_data_step2, 
                                                        all_pseudogene_zone_3_OR_genes_pie_data_step2, 
                                                        all_pseudogene_zone_4_OR_genes_pie_data_step2, 
                                                        all_pseudogene_zone_5_OR_genes_pie_data_step2),
                                      changes=rep(c("Upregulated", "Downregulated", "Not_changed"), 18),
                                      gene_features=c(rep(paste("all ", all_genes_number_step2, " genes", sep = ""), 3), 
                                                      rep(paste(all_functional_TAAR_genes_number_step2, " functional TAAR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_OR_genes_number_step2, " functional OR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_ClassI_OR_genes_number_step2, " functional Class I OR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_ClassII_OR_genes_number_step2, " functional Class II OR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_zone_1_OR_genes_number_step2, " functional zone 1 OR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_zone_2_OR_genes_number_step2, " functional zone 2 OR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_zone_3_OR_genes_number_step2, " functional zone 3 OR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_zone_4_OR_genes_number_step2, " functional zone 4 OR genes", sep = ""), 3), 
                                                      rep(paste(all_functional_zone_5_OR_genes_number_step2, " functional zone 5 OR genes", sep = ""), 3), 
                                                      rep(paste(all_pseudogene_OR_genes_number_step2, " OR pseudogenes", sep = ""), 3), 
                                                      rep(paste(all_pseudogene_ClassI_OR_genes_number_step2, " Class I OR pseudogenes", sep = ""), 3), 
                                                      rep(paste(all_pseudogene_ClassII_OR_genes_number_step2, " Class II OR pseudogenes", sep = ""), 3),
                                                      rep(paste(all_pseudogene_zone_1_OR_genes_number_step2, " zone 1 pseudogene OR genes", sep = ""), 3), 
                                                      rep(paste(all_pseudogene_zone_2_OR_genes_number_step2, " zone 2 pseudogene OR genes", sep = ""), 3), 
                                                      rep(paste(all_pseudogene_zone_3_OR_genes_number_step2, " zone 3 pseudogene OR genes", sep = ""), 3), 
                                                      rep(paste(all_pseudogene_zone_4_OR_genes_number_step2, " zone 4 pseudogene OR genes", sep = ""), 3), 
                                                      rep(paste(all_pseudogene_zone_5_OR_genes_number_step2, " zone 5 pseudogene OR genes", sep = ""), 3) 
                                      ))
all_data_bar_plot_step2$changes <- factor(all_data_bar_plot_step2$changes, levels=c("Downregulated", "Upregulated", "Not_changed")) #change to factor so that we can change the stack order

plot_step2_1 <- ggplot(all_data_bar_plot_step2, aes(x = gene_features, y = Number_of_genes, fill = changes)) + 
  geom_col(position = "fill") + 
  scale_fill_manual(values = c("springgreen3", "lightcoral", "gray80")) + #change the color
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(limits = c(paste("all ", all_genes_number_step2, " genes", sep = ""),
                              paste(all_functional_TAAR_genes_number_step2, " functional TAAR genes", sep = ""),
                              paste(all_functional_OR_genes_number_step2, " functional OR genes", sep = ""),
                              paste(all_functional_ClassI_OR_genes_number_step2, " functional Class I OR genes", sep = ""),
                              paste(all_functional_ClassII_OR_genes_number_step2, " functional Class II OR genes", sep = ""),
                              paste(all_functional_zone_1_OR_genes_number_step2, " functional zone 1 OR genes", sep = ""),
                              paste(all_functional_zone_2_OR_genes_number_step2, " functional zone 2 OR genes", sep = ""),
                              paste(all_functional_zone_3_OR_genes_number_step2, " functional zone 3 OR genes", sep = ""),
                              paste(all_functional_zone_4_OR_genes_number_step2, " functional zone 4 OR genes", sep = ""),
                              paste(all_functional_zone_5_OR_genes_number_step2, " functional zone 5 OR genes", sep = ""),
                              paste(all_pseudogene_OR_genes_number_step2, " OR pseudogenes", sep = ""), 
                              paste(all_pseudogene_ClassI_OR_genes_number_step2, " Class I OR pseudogenes", sep = ""),
                              paste(all_pseudogene_ClassII_OR_genes_number_step2, " Class II OR pseudogenes", sep = ""),
                              paste(all_pseudogene_zone_1_OR_genes_number_step2, " zone 1 pseudogene OR genes", sep = ""),
                              paste(all_pseudogene_zone_2_OR_genes_number_step2, " zone 2 pseudogene OR genes", sep = ""),
                              paste(all_pseudogene_zone_3_OR_genes_number_step2, " zone 3 pseudogene OR genes", sep = ""),
                              paste(all_pseudogene_zone_4_OR_genes_number_step2, " zone 4 pseudogene OR genes", sep = ""),
                              paste(all_pseudogene_zone_5_OR_genes_number_step2, " zone 5 pseudogene OR genes", sep = ""))) + # change the order of x axis
  ggtitle('Trim66_KO homo vs het with threshold padj < 0.05 and abs(FoldChange)>= 1.5 to determine downregulated and upregulated genes') +
  theme(plot.title = element_text(hjust=0.3, vjust = 1, size = 10)) + # adjust the title position and size
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

pdf(file = paste(resultDir,"/stack_bar_plots_changed_genes_percentage_Trim66_KO_homo_vs_het.pdf", sep=''), width=30/2.54, height = 15/2.54) # width and height in inches
plot_step2_1
ggplot() +  #check https://statisticsglobe.com/plot-only-text-in-r
  annotate("text",
           x = 0,
           y = 0,
           size = 1.5, #adjust the size if the data are not fully printed
           label = paste(capture.output(all_data_bar_plot_step2), collapse = "\n")) + # check https://stackoverflow.com/questions/14326573/print-a-data-frame-in-the-white-space-of-a-plot
  theme_void()
dev.off()


#####################################################################################################################
# save the results of OR, TAAR, and Ms4a
res_OR <- res_all[grepl('^Olfr', rownames(res_all)), ]
write.csv(res_OR, file = paste(baseDir, "/results/results_OR_genes_Gencode_Trim66_homo_vs_het.csv", sep=''))
res_TAAR <- res_all[grepl('^Taar', rownames(res_all)), ]
write.csv(res_TAAR, file = paste(baseDir, "/results/results_TAAR_genes_Gencode_Trim66_homo_vs_het.csv", sep=''))
res_Ms4a <- res_all[grepl('^Ms4a', rownames(res_all)), ]
write.csv(res_Ms4a, file = paste(baseDir, "/results/results_Ms4a_genes_Gencode_Trim66_homo_vs_het.csv", sep=''))
# functional OR genes are defined as: 1. having “Olfr” in the gene name; 2. having “protein_coding” in the gene type
# so there might be "ps" in gene name but actually are functional OR genes. check the paper: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6583-3#MOESM2
res_functional_OR <- res_all[grepl('^Olfr', rownames(res_all)) & grepl('protein_coding', res_all$gene_type), ]
write.csv(res_functional_OR, file = paste(baseDir, "/results/results_functional_OR_genes_Gencode_Trim66_homo_vs_het.csv", sep=''))
res_pseudogene_OR <- res_all[grepl('^Olfr', rownames(res_all)) & grepl('pseudogene', res_all$gene_type), ]
write.csv(res_pseudogene_OR, file = paste(baseDir, "/results/results_pseudogene_OR_genes_Gencode_Trim66_homo_vs_het.csv", sep=''))
dim(res_OR)
dim(res_functional_OR)
dim(res_pseudogene_OR)

# Plot step 3: to check if ORs and TAARs in different cluster and if the distance to OR/TAAR enhancers matter ---------------------------

#####################################################################################################################
     # plot step 3: to check if ORs and TAARs in different cluster and if the distance to OR/TAAR enhancers matter
#####################################################################################################################
#3.1 functional ORs and TAARs along real chromosome coordinates
###################################################################################################################### add one column to indicate if this receptor is significantly chaged in homo vs het
res_functional_OR_and_TAAR_step3.1 <- res_all[grepl("Olfr|Taar[2-9]", rownames(res_all)) & grepl("protein_coding", res_all$gene_type), ]
res_functional_OR_and_TAAR_step3.1$if_significant_homoVShet <- res_functional_OR_and_TAAR_step3.1$padj<0.05 & abs(res_functional_OR_and_TAAR_step3.1$log2FoldChange)>=0.585
# replace NA with FALSE
res_functional_OR_and_TAAR_step3.1[is.na(res_functional_OR_and_TAAR_step3.1$padj), ]$if_significant_homoVShet <- FALSE
# add one column to indicate if this receptor is significantly upregulated or downregulated in homo vs het
# note that the code in step3.5 is better. Use that code next time. This is fine for small number of genes and is for record.
res_functional_OR_and_TAAR_step3.1$if_significant_upregulate_or_downregulate <- NA
for (i in 1:nrow(res_functional_OR_and_TAAR_step3.1)) {
  if (is.na(res_functional_OR_and_TAAR_step3.1[i, ]$padj)) {
    res_functional_OR_and_TAAR_step3.1[i, ]$if_significant_upregulate_or_downregulate <- "Not_changed"
  } else if (res_functional_OR_and_TAAR_step3.1[i, ]$padj <0.05 & res_functional_OR_and_TAAR_step3.1[i, ]$log2FoldChange >=0.585) {
    res_functional_OR_and_TAAR_step3.1[i, ]$if_significant_upregulate_or_downregulate <- "Upregulated_significantly"
  } else if (res_functional_OR_and_TAAR_step3.1[i, ]$padj <0.05 & res_functional_OR_and_TAAR_step3.1[i, ]$log2FoldChange <=-0.585) {
    res_functional_OR_and_TAAR_step3.1[i, ]$if_significant_upregulate_or_downregulate <- "Downregulated_significantly"
  } else {
    res_functional_OR_and_TAAR_step3.1[i, ]$if_significant_upregulate_or_downregulate <- "Not_changed"
  }
}
# add one column to inverse the chromosome numbers to negative numbers
unique(res_functional_OR_and_TAAR_step3.1$chr_name)
res_functional_OR_and_TAAR_step3.1$reversed_chr_number <- ifelse(res_functional_OR_and_TAAR_step3.1$chr_name == "chrX", -20, as.numeric(substring(res_functional_OR_and_TAAR_step3.1$chr_name, 4)) * -1) #check substring function, https://stackoverflow.com/questions/17215789/extract-a-substring-according-to-a-pattern
# plot the ORs based on chromosome numbers, and strand direction
res_functional_OR_and_TAAR_dataframe_step3.1 <- as.data.frame(res_functional_OR_and_TAAR_step3.1)
plot_step3.1 <- ggplot(data=res_functional_OR_and_TAAR_dataframe_step3.1, aes(x=genomic_start, y=reversed_chr_number, 
                                                                              colour=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate,
                                                                              alpha=if_significant_upregulate_or_downregulate, shape=as.factor(genomic_strand))) +
  geom_point(size=2.5, position = position_jitter(width = 0, height = 0.4)) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_shape_manual(values = c("\u25C4","\u25BA")) + #to choose other shapes like different triangles. https://stackoverflow.com/questions/30742379/creating-new-shape-palettes-in-ggplot2-and-other-r-graphics
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,1,1)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("chr coordinate") + ylab("chr number") +
  xlim(0, 2e8) + # 2.0e+8 as the maximum to compare plot_step3.1 with the following plot_step3.5
  scale_y_continuous(breaks = -1:-20, labels = c(1:19, "X"))  + # set the breaks to show as indicated, and also change the labels as indicated
  labs(colour = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "\u25BA plus strand\n\u25C4 minus strand") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
# The geometric shape may look weird if you print to pdf file. check this website for solution:
# https://stackoverflow.com/questions/25053072/how-to-draw-half-filled-points-in-r-preferably-using-ggplot/34083063#34083063
grDevices::cairo_pdf(paste(resultDir, "/functional_ORs_TAARs_changed_in_genome_Trim66_KO_whole_genome.pdf", sep=""), width = 21, height = 27.9, family="Arial Unicode MS")
plot_step3.1
dev.off()

# step3.1_2: add the bar of 63 OR enhancers, 1 J element, and 2 TAAR enhancers
enhancers_info <- read.delim("/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/OR_63_enhancers_and_1_J_element_and_2_TAAR_enhancers.bed", header=F, sep="\t",
                             row.names=NULL, as.is=T)
colnames(enhancers_info) <- c("chr_name", "genomic_start", "genomic_end", "enhancer_name")
enhancers_info$reversed_chr_number <- ifelse(enhancers_info$chr_name == "chrX", -20, as.numeric(substring(enhancers_info$chr_name, 4)) * -1)
plot_step3.1_2_with_enhancer_info <- ggplot() +
  geom_segment(data = enhancers_info,
               aes(x=genomic_start, y=reversed_chr_number - 0.4, xend=genomic_start, yend=reversed_chr_number + 0.4), 
               color="gray60", 
               alpha=1,
               size=0.1)+
  geom_text(data = enhancers_info,  
            aes(x=genomic_start,
                y=reversed_chr_number + 0.5,
                label = enhancer_name), 
            size =1, colour='black', 
            hjust=0, #left aligned
            position=position_jitter(height=0.1)) +
  geom_point(data=res_functional_OR_and_TAAR_dataframe_step3.1, aes(x=genomic_start, y=reversed_chr_number, 
                                                                    colour=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate,
                                                                    alpha=if_significant_upregulate_or_downregulate, shape=as.factor(genomic_strand)), 
             size=0.25, position = position_jitter(width = 0, height = 0.4)) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_shape_manual(values = c("\u25C4","\u25BA")) + #to choose other shapes like different triangles. https://stackoverflow.com/questions/30742379/creating-new-shape-palettes-in-ggplot2-and-other-r-graphics
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,1,1)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 \nwith 63 OR enhancers, 1 J elements, and 2 TAAR enhancers as grey bars') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("chr coordinate") + ylab("chr number") +
  xlim(0, 2e8) + # 2.0e+8 as the maximum to compare plot_step3.1 with the following plot_step3.5
  scale_y_continuous(breaks = -1:-20, labels = c(1:19, "X"))  + # set the breaks to show as indicated, and also change the labels as indicated
  labs(colour = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "\u25BA plus strand\n\u25C4 minus strand") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
# save the plot as pdf
grDevices::cairo_pdf(paste(resultDir, "/functional_ORs_TAARs_changed_in_genome_Trim66_KO_whole_genome_with_enhancers.pdf", sep=""), width = 21, height = 27.9, family="Arial Unicode MS")
plot_step3.1_2_with_enhancer_info
dev.off()

# step3.1_3: plot the log2FC against the distance of functional ORs and TAARs to the closest enhancers
#####################################################################################################################
# Plot Step 3.1_3: plot OR and TAAR changes with their distance to enhancers to see if they are correlated
#####################################################################################################################
#3.1_3 functional ORs and TAARs
res_functional_OR_and_TAAR_dataframe_step3.1_3 <- res_functional_OR_and_TAAR_dataframe_step3.1
# add one column to indicate TAAR or OR
res_functional_OR_and_TAAR_dataframe_step3.1_3$TAAR_or_OR <- ifelse(grepl("^Taar", rownames(res_functional_OR_and_TAAR_dataframe_step3.1_3)), "TAAR", "OR")
# add one column to change NA in log2FoldChange to 0 so that those data will be not removed when plotting
res_functional_OR_and_TAAR_dataframe_step3.1_3$log2FoldChange_NA_to_0 <- res_functional_OR_and_TAAR_dataframe_step3.1_3$log2FoldChange
res_functional_OR_and_TAAR_dataframe_step3.1_3[is.na(res_functional_OR_and_TAAR_dataframe_step3.1_3$log2FoldChange),]$log2FoldChange_NA_to_0 <- 0
# add one column to show the distances to the closest enhancers. The distances are between the mean(genomic_start, genomic_end) of ORs and the mean(genomic_start, genomic_end) of enhancers.
res_functional_OR_and_TAAR_dataframe_step3.1_3$distance_to_closest_enhancer <- 0
res_functional_OR_and_TAAR_dataframe_step3.1_3$the_closest_enhancer_name <- ""
for (i in 1:nrow(res_functional_OR_and_TAAR_dataframe_step3.1_3)) {
  index_3.1_3 <- which(enhancers_info$chr_name == res_functional_OR_and_TAAR_dataframe_step3.1_3[i, "chr_name"])
  res_functional_OR_and_TAAR_dataframe_step3.1_3[i,]$distance_to_closest_enhancer <- 
    min(abs((res_functional_OR_and_TAAR_dataframe_step3.1_3[i, "genomic_start"] + res_functional_OR_and_TAAR_dataframe_step3.1_3[i,"genomic_end"])/2 -
              (enhancers_info[index_3.1_3, "genomic_start"] + enhancers_info[index_3.1_3, "genomic_end"])/2))
  if (length(index_3.1_3) != 0) {
    enhancers_temp_3.1_3 <- enhancers_info[index_3.1_3, ]
    res_functional_OR_and_TAAR_dataframe_step3.1_3[i,]$the_closest_enhancer_name <- 
      enhancers_temp_3.1_3[which.min(abs((res_functional_OR_and_TAAR_dataframe_step3.1_3[i, "genomic_start"] + res_functional_OR_and_TAAR_dataframe_step3.1_3[i,"genomic_end"])/2 -
                                           (enhancers_temp_3.1_3$genomic_start + enhancers_temp_3.1_3$genomic_end)/2)), ]$enhancer_name
  }
}
# There are warning message that "returning Inf" because there are no enhancers in chr5 and chr8 for functional ORs
# What are they? It turns out that there are 1 OR in chr5 and 5 ORs in chr8 showing Inf in distance_to_closest_enhancer
res_functional_OR_and_TAAR_dataframe_step3.1_3[which(is.infinite(res_functional_OR_and_TAAR_dataframe_step3.1_3$distance_to_closest_enhancer)), ]
# delete the 6 ORs in chr5 and chr8
dim(res_functional_OR_and_TAAR_dataframe_step3.1_3)
res_functional_OR_and_TAAR_dataframe_step3.1_3 <- res_functional_OR_and_TAAR_dataframe_step3.1_3[is.finite(res_functional_OR_and_TAAR_dataframe_step3.1_3$distance_to_closest_enhancer), ]
dim(res_functional_OR_and_TAAR_dataframe_step3.1_3)

# plot the receptor changes against distance to the closest enhancers
plot_step3.1_3_a <- ggplot(data=res_functional_OR_and_TAAR_dataframe_step3.1_3, aes(x=distance_to_closest_enhancer, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR), #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR),
                 alpha=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR)), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,1,0.8,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 10000000, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_functional_OR_and_TAAR_dataframe_step3.1_3$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step3.1_3$distance_to_closest_enhancer), 4),
                         ", p value: ",
                         round(cor.test(res_functional_OR_and_TAAR_dataframe_step3.1_3$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step3.1_3$distance_to_closest_enhancer)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Distance between ORs/TAARs and the closest enhancers (middle of ORs/TAARs to middle of enhancer)") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
res_functional_OR_and_TAAR_dataframe_step3.1_3$name <- rownames(res_functional_OR_and_TAAR_dataframe_step3.1_3)
plot_step3.1_3_b <- ggplot(data=res_functional_OR_and_TAAR_dataframe_step3.1_3, aes(x=distance_to_closest_enhancer, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR), #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR),
                 alpha=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR)), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_alpha_manual(values = c(0.5,0.5,0.5,1,0.8,1)) +
  geom_text(data=subset(res_functional_OR_and_TAAR_dataframe_step3.1_3, (padj<0.05 & abs(log2FoldChange)>=0.585) ==TRUE), 
            aes(label = name), 
            size =5, colour='blue') +  # label=name in the aes() suggest the text to label
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 10000000, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_functional_OR_and_TAAR_dataframe_step3.1_3$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step3.1_3$distance_to_closest_enhancer), 4),
                         ", p value: ",
                         round(cor.test(res_functional_OR_and_TAAR_dataframe_step3.1_3$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step3.1_3$distance_to_closest_enhancer)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Distance between ORs/TAARs and the closest enhancers (middle of ORs/TAARs to middle of enhancer)") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

arrange_step3.1_3 <- ggarrange(plot_step3.1_3_a, plot_step3.1_3_b, ncol = 1, nrow = 2)
ggsave(paste(baseDir, "/results/functional_ORs_TAARs_changes_against_distance_to_enhancers_Trim66_KO.pdf", sep=''), arrange_step3.1_3, width = 120, height = 60, units ="cm")

# plot the mean+-se of distance to the closest enhancers. The distances are between the mean(genomic_start, genomic_end) of ORs and the mean(genomic_start, genomic_end) of enhancers.
# obtain the summary data
res_functional_OR_and_TAAR_dataframe_grouped_summary_step3.1_3 <- res_functional_OR_and_TAAR_dataframe_step3.1_3 %>% 
  mutate(if_significant_upregulate_or_downregulate = factor(if_significant_upregulate_or_downregulate, 
                                                            levels = c("Downregulated_significantly", "Upregulated_significantly", "Not_changed"))) %>% 
  group_by(if_significant_upregulate_or_downregulate) %>% 
  summarize(mean_distance_to_closest_enhancer=mean(distance_to_closest_enhancer), 
            sd_distance_to_closest_enhancer=sd(distance_to_closest_enhancer), 
            Numbers=n(), 
            se=sd_distance_to_closest_enhancer/sqrt(mean_distance_to_closest_enhancer), 
            upper_limit=mean_distance_to_closest_enhancer+se, 
            lower_limit=mean_distance_to_closest_enhancer-se)

# obtain the one way ANOVA test results
one_way_ANOVA_for_OR_and_TAAR_step3.1_3 <- aov(distance_to_closest_enhancer ~ if_significant_upregulate_or_downregulate, data = res_functional_OR_and_TAAR_dataframe_step3.1_3) # one way ANOVA as there are 3 types, check: https://www.scribbr.com/statistics/anova-in-r/
summary(one_way_ANOVA_for_OR_and_TAAR_step3.1_3)
tukey_post_hoc_for_OR_and_TAAR_step3.1_3<-TukeyHSD(one_way_ANOVA_for_OR_and_TAAR_step3.1_3)
tukey_post_hoc_for_OR_and_TAAR_step3.1_3

# start plotting
plot_step3.1_3_c <- ggplot() + 
  geom_bar(res_functional_OR_and_TAAR_dataframe_grouped_summary_step3.1_3,
           mapping=aes(x=if_significant_upregulate_or_downregulate, y=mean_distance_to_closest_enhancer, color=if_significant_upregulate_or_downregulate),
           stat="identity", 
           position = position_dodge(0.95),
           width = 0.9, #adjust the bar width
           size=0.25, # adjust the line thickness
           fill="gray95") + 
  geom_errorbar(res_functional_OR_and_TAAR_dataframe_grouped_summary_step3.1_3,
                mapping=aes(x=if_significant_upregulate_or_downregulate, ymin=lower_limit, ymax=upper_limit, color=if_significant_upregulate_or_downregulate),
                position = position_dodge(0.95), width = .2, size=0.25) +
  scale_color_manual(values = c("springgreen3", "red1", "black")) +
  geom_point(res_functional_OR_and_TAAR_dataframe_step3.1_3,
             mapping=aes(x=if_significant_upregulate_or_downregulate, 
                         y=distance_to_closest_enhancer, 
                         group=if_significant_upregulate_or_downregulate, 
                         fill=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR)),
             position=position_jitterdodge(jitter.width = 3,
                                           dodge.width = 3),
             size=3, shape=21, stroke = 0.4, # stroke control border thickness
             color="black", alpha =0.5) +
  scale_fill_manual(values = c(rep("NA", 3), "blue1", "black", "magenta1")) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(2)), # adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5))) +
  ggtitle('Distance between functional ORs/TAARs and the closest enhancers (middle of ORs/TAARs to middle of enhancer)') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("significantly changed?") + 
  ylab("distance between ORs/TAARs and the closest enhancer (mean +- se)") +
  ggbreak::scale_y_break(c(5000000, 36000000)) +# check https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html#introduction
  annotate("text", x = 2, y = 4000000, 
           label = paste(capture.output(res_functional_OR_and_TAAR_dataframe_grouped_summary_step3.1_3), sep=" ",collapse = "\n"),
           size=4, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  annotate("text", x = 2, y = 3000000, 
           label = "One way ANOVA and post hoc Tukey analysis",
           size=4, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  annotate("text", x = 2, y = 2000000, 
           label = paste(capture.output(tukey_post_hoc_for_OR_and_TAAR_step3.1_3), sep=" ",collapse = "\n"),
           size=4, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

plot_step3.1_3_c
ggsave(paste(baseDir, "/results/functional_ORs_TAARs_distance_to_enhancers_Trim66_KO.pdf", sep=''), plot_step3.1_3_c, width = 50, height = 25, units ="cm")


#####################################################################################################################
#3.2 pseudogene ORs along real chromosome coordinates
###################################################################################################################### add one column to indicate if this receptor is significantly chaged in homo vs het
res_pseudogene_OR_step3.2 <- res_pseudogene_OR
res_pseudogene_OR_step3.2$if_significant_homoVShet <- res_pseudogene_OR_step3.2$padj<0.05 & abs(res_pseudogene_OR_step3.2$log2FoldChange)>=0.585
# replace NA with FALSE
res_pseudogene_OR_step3.2[is.na(res_pseudogene_OR_step3.2$padj), ]$if_significant_homoVShet <- FALSE
# add one column to indicate if this receptor is significantly upregulated or downregulated in homo vs het
# note that the code in step3.5 is better. Use that code next time. This is fine for small number of genes and is for record.
res_pseudogene_OR_step3.2$if_significant_upregulate_or_downregulate <- NA
for (i in 1:nrow(res_pseudogene_OR_step3.2)) {
  if (is.na(res_pseudogene_OR_step3.2[i, ]$padj)) {
    res_pseudogene_OR_step3.2[i, ]$if_significant_upregulate_or_downregulate <- "Not_changed"
  } else if (res_pseudogene_OR_step3.2[i, ]$padj <0.05 & res_pseudogene_OR_step3.2[i, ]$log2FoldChange >=0.585) {
    res_pseudogene_OR_step3.2[i, ]$if_significant_upregulate_or_downregulate <- "Upregulated_significantly"
  } else if (res_pseudogene_OR_step3.2[i, ]$padj <0.05 & res_pseudogene_OR_step3.2[i, ]$log2FoldChange <=-0.585) {
    res_pseudogene_OR_step3.2[i, ]$if_significant_upregulate_or_downregulate <- "Downregulated_significantly"
  } else {
    res_pseudogene_OR_step3.2[i, ]$if_significant_upregulate_or_downregulate <- "Not_changed"
  }
}
# add one column to inverse the chromosome numbers to negative numbers
res_pseudogene_OR_step3.2$reversed_chr_number <- ifelse(res_pseudogene_OR_step3.2$chr_name == "chrX", -20, as.numeric(substring(res_pseudogene_OR_step3.2$chr_name, 4)) * -1) #check substring function, https://stackoverflow.com/questions/17215789/extract-a-substring-according-to-a-pattern
# plot the ORs based on chromosome numbers, and strand direction
res_pseudogene_OR_dataframe_step3.2 <- as.data.frame(res_pseudogene_OR_step3.2)
plot_step3.2 <- ggplot(data=res_pseudogene_OR_dataframe_step3.2, aes(x=genomic_start, y=reversed_chr_number, 
                                                                     colour=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate,
                                                                     alpha=if_significant_upregulate_or_downregulate, shape=as.factor(genomic_strand))) +
  geom_point(size=2.5, position = position_jitter(width = 0, height = 0.4)) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_shape_manual(values = c("\u25C4","\u25BA")) + #to choose other shapes like different triangles. https://stackoverflow.com/questions/30742379/creating-new-shape-palettes-in-ggplot2-and-other-r-graphics
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,1,1)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)),# adjust the title size
        axis.title = element_text(size = rel(1.5)), 
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("chr coordinate") + ylab("chr number") +
  xlim(0, 2e8) + # 2.0e+8 as the maximum to compare plot_step3.2 with the following plot_step3.5
  scale_y_continuous(breaks = -1:-20, labels = c(1:19, "X"))  + # set the breaks to show as indicated, and also change the labels as indicated
  labs(colour = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "\u25BA plus strand\n\u25C4 minus strand") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
# The geometric shape may look weird if you print to pdf file. check this website for solution:
# https://stackoverflow.com/questions/25053072/how-to-draw-half-filled-points-in-r-preferably-using-ggplot/34083063#34083063
grDevices::cairo_pdf(paste(resultDir, "/pseudogenes_ORs_changed_in_genome_Trim66_KO_whole_genome.pdf", sep=""), width = 21, height = 27.9, family="Arial Unicode MS")
plot_step3.2
dev.off()

# add the bar of 63 OR enhancers, 1 J element, and 2 TAAR enhancers
plot_step3.2_2_with_enhancer_info <- ggplot() +
  geom_segment(data = enhancers_info,
               aes(x=genomic_start, y=reversed_chr_number - 0.4, xend=genomic_start, yend=reversed_chr_number + 0.4), 
               color="gray60", 
               alpha=1,
               size=0.1)+
  geom_text(data = enhancers_info,  
            aes(x=genomic_start,
                y=reversed_chr_number + 0.5,
                label = enhancer_name), 
            size =1, colour='black', 
            hjust=0, #left aligned
            position=position_jitter(height=0.1)) +
  geom_point(data=res_pseudogene_OR_dataframe_step3.2, aes(x=genomic_start, y=reversed_chr_number, 
                                                                    colour=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate,
                                                                    alpha=if_significant_upregulate_or_downregulate, shape=as.factor(genomic_strand)), 
             size=0.25, position = position_jitter(width = 0, height = 0.4)) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_shape_manual(values = c("\u25C4","\u25BA")) + #to choose other shapes like different triangles. https://stackoverflow.com/questions/30742379/creating-new-shape-palettes-in-ggplot2-and-other-r-graphics
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,1,1)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 \nwith 63 OR enhancers, 1 J elements, and 2 TAAR enhancers as grey bars') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("chr coordinate") + ylab("chr number") +
  xlim(0, 2e8) + # 2.0e+8 as the maximum to compare plot_step3.1 with the following plot_step3.5
  scale_y_continuous(breaks = -1:-20, labels = c(1:19, "X"))  + # set the breaks to show as indicated, and also change the labels as indicated
  labs(colour = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "\u25BA plus strand\n\u25C4 minus strand") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
# save the plot as pdf
grDevices::cairo_pdf(paste(resultDir, "/pseudogenes_ORs_changed_in_genome_Trim66_KO_whole_genome_with_enhancers.pdf", sep=""), width = 21, height = 27.9, family="Arial Unicode MS")
plot_step3.2_2_with_enhancer_info
dev.off()

# step3.2_3: plot the log2FC against the distance of pseudogene ORs to the closest enhancers
#####################################################################################################################
# Plot Step 3.2_3: plot pseudogene ORs changes with their distance to enhancers to see if they are correlated
#####################################################################################################################
#3.2_3 pseudogene ORs
res_pseudogene_OR_dataframe_step3.2_3 <- res_pseudogene_OR_dataframe_step3.2
# add one column to change NA in log2FoldChange to 0 so that those data will be not removed when plotting
res_pseudogene_OR_dataframe_step3.2_3$log2FoldChange_NA_to_0 <- res_pseudogene_OR_dataframe_step3.2_3$log2FoldChange
res_pseudogene_OR_dataframe_step3.2_3[is.na(res_pseudogene_OR_dataframe_step3.2_3$log2FoldChange),]$log2FoldChange_NA_to_0 <- 0
# add one column to show the distances to the closest enhancers. The distances are between the mean(genomic_start, genomic_end) of ORs and the mean(genomic_start, genomic_end) of enhancers.
res_pseudogene_OR_dataframe_step3.2_3$distance_to_closest_enhancer <- 0
res_pseudogene_OR_dataframe_step3.2_3$the_closest_enhancer_name <- ""
for (i in 1:nrow(res_pseudogene_OR_dataframe_step3.2_3)) {
  index_3.2_3 <- which(enhancers_info$chr_name == res_pseudogene_OR_dataframe_step3.2_3[i, "chr_name"])
  res_pseudogene_OR_dataframe_step3.2_3[i,]$distance_to_closest_enhancer <- 
    min(abs((res_pseudogene_OR_dataframe_step3.2_3[i, "genomic_start"] + res_pseudogene_OR_dataframe_step3.2_3[i,"genomic_end"])/2 -
              (enhancers_info[index_3.2_3, "genomic_start"] + enhancers_info[index_3.2_3, "genomic_end"])/2))
  if (length(index_3.2_3) != 0) {
    enhancers_temp_3.2_3 <- enhancers_info[index_3.2_3, ]
    res_pseudogene_OR_dataframe_step3.2_3[i,]$the_closest_enhancer_name <- 
      enhancers_temp_3.2_3[which.min(abs((res_pseudogene_OR_dataframe_step3.2_3[i, "genomic_start"] + res_pseudogene_OR_dataframe_step3.2_3[i,"genomic_end"])/2 -
                                           (enhancers_temp_3.2_3$genomic_start + enhancers_temp_3.2_3$genomic_end)/2)), ]$enhancer_name
  }
}
# There are warning message that "returning Inf" because there are no enhancers in chr5 and chr8 for functional ORs
# What are they? It turns out that there are 1 OR in chr5 and 5 ORs in chr8 showing Inf in distance_to_closest_enhancer
res_pseudogene_OR_dataframe_step3.2_3[which(is.infinite(res_pseudogene_OR_dataframe_step3.2_3$distance_to_closest_enhancer)), ]
# delete the 6 ORs in chr5 and chr8
dim(res_pseudogene_OR_dataframe_step3.2_3)
res_pseudogene_OR_dataframe_step3.2_3 <- res_pseudogene_OR_dataframe_step3.2_3[is.finite(res_pseudogene_OR_dataframe_step3.2_3$distance_to_closest_enhancer), ]
dim(res_pseudogene_OR_dataframe_step3.2_3)

# plot the receptor changes against distance to the closest enhancers
plot_step3.2_3_a <- ggplot(data=res_pseudogene_OR_dataframe_step3.2_3, aes(x=distance_to_closest_enhancer, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate,
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1,1,0.8,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 30000000, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_pseudogene_OR_dataframe_step3.2_3$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step3.2_3$distance_to_closest_enhancer), 4),
                         ", p value: ",
                         round(cor.test(res_pseudogene_OR_dataframe_step3.2_3$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step3.2_3$distance_to_closest_enhancer)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Distance between pseudogene ORs and the closest enhancers (middle of pesudogene ORs to middle of enhancer)") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
res_pseudogene_OR_dataframe_step3.2_3$name <- rownames(res_pseudogene_OR_dataframe_step3.2_3)
plot_step3.2_3_b <- ggplot(data=res_pseudogene_OR_dataframe_step3.2_3, aes(x=distance_to_closest_enhancer, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate,
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1,1,0.8,1)) +
  geom_text(data=subset(res_pseudogene_OR_dataframe_step3.2_3, (padj<0.05 & abs(log2FoldChange)>=0.585) ==TRUE), 
            aes(label = name), 
            size =5, colour='blue') +  # label=name in the aes() suggest the text to label
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 30000000, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_pseudogene_OR_dataframe_step3.2_3$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step3.2_3$distance_to_closest_enhancer), 4),
                         ", p value: ",
                         round(cor.test(res_pseudogene_OR_dataframe_step3.2_3$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step3.2_3$distance_to_closest_enhancer)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Distance between pseudogene ORs and the closest enhancers (middle of pesudogene ORs to middle of enhancer)") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

arrange_step3.2_3 <- ggarrange(plot_step3.2_3_a, plot_step3.2_3_b, ncol = 1, nrow = 2)
ggsave(paste(baseDir, "/results/pseudogene_ORs_changes_against_distance_to_enhancers_Trim66_KO.pdf", sep=''), arrange_step3.2_3, width = 120, height = 60, units ="cm")

# plot the mean of distance to the closest enhancers. The distances are between the mean(genomic_start, genomic_end) of ORs and the mean(genomic_start, genomic_end) of enhancers.
# obtain the summary data
res_pseudogene_OR_dataframe_grouped_summary_step3.2_3 <- res_pseudogene_OR_dataframe_step3.2_3 %>% 
  mutate(if_significant_upregulate_or_downregulate = factor(if_significant_upregulate_or_downregulate, 
                                                            levels = c("Downregulated_significantly", "Upregulated_significantly", "Not_changed"))) %>% 
  group_by(if_significant_upregulate_or_downregulate) %>% 
  summarize(mean_distance_to_closest_enhancer=mean(distance_to_closest_enhancer), 
            sd_distance_to_closest_enhancer=sd(distance_to_closest_enhancer), 
            Numbers=n(), 
            se=sd_distance_to_closest_enhancer/sqrt(mean_distance_to_closest_enhancer), 
            upper_limit=mean_distance_to_closest_enhancer+se, 
            lower_limit=mean_distance_to_closest_enhancer-se)

# obtain the one way ANOVA test results
one_way_ANOVA_for_pseudogene_OR_step3.2_3 <- aov(distance_to_closest_enhancer ~ if_significant_upregulate_or_downregulate, data = res_pseudogene_OR_dataframe_step3.2_3) # one way ANOVA as there are 3 types, check: https://www.scribbr.com/statistics/anova-in-r/
summary(one_way_ANOVA_for_pseudogene_OR_step3.2_3)
tukey_post_hoc_for_pseudogene_OR_step3.2_3<-TukeyHSD(one_way_ANOVA_for_pseudogene_OR_step3.2_3)
tukey_post_hoc_for_pseudogene_OR_step3.2_3

# start plotting
plot_step3.2_3_c <- ggplot() + 
  geom_bar(res_pseudogene_OR_dataframe_grouped_summary_step3.2_3,
           mapping=aes(x=if_significant_upregulate_or_downregulate, y=mean_distance_to_closest_enhancer, color=if_significant_upregulate_or_downregulate),
           stat="identity", 
           position = position_dodge(0.95),
           width = 0.9, #adjust the bar width
           size=0.25, # adjust the line thickness
           fill="gray95") + 
  geom_errorbar(res_pseudogene_OR_dataframe_grouped_summary_step3.2_3,
                mapping=aes(x=if_significant_upregulate_or_downregulate, ymin=lower_limit, ymax=upper_limit, color=if_significant_upregulate_or_downregulate),
                position = position_dodge(0.95), width = .2, size=0.25) +
  scale_color_manual(values = c("springgreen3", "red1", "black")) +
  geom_point(res_pseudogene_OR_dataframe_step3.2_3,
             mapping=aes(x=if_significant_upregulate_or_downregulate, y=distance_to_closest_enhancer, group=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate),
             position=position_jitterdodge(jitter.width = 0.75,
                                           dodge.width = 0.9),
             size=3, shape=21, stroke = 0.4, # stroke control border thickness
             color="black", alpha =0.5) +
  scale_fill_manual(values = rep("white", 3)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(2)), # adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5))) +
  ggtitle('Distance between pseudogene ORs and the closest enhancers (middle of pesudogene ORs to middle of enhancer)') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("significantly changed?") + 
  ylab("distance between pseudogene ORs and the closest enhancer (mean +- se)") +
  ggbreak::scale_y_break(c(7500000, 90000000)) +# check https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html#introduction
  annotate("text", x = 2, y = 7000000, 
           label = paste(capture.output(res_pseudogene_OR_dataframe_grouped_summary_step3.2_3), sep=" ",collapse = "\n"),
           size=4, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  annotate("text", x = 2, y = 5000000, 
           label = "One way ANOVA and post hoc Tukey analysis",
           size=4, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  annotate("text", x = 2, y = 3000000, 
           label = paste(capture.output(tukey_post_hoc_for_pseudogene_OR_step3.2_3), sep=" ",collapse = "\n"),
           size=4, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

plot_step3.2_3_c
ggsave(paste(baseDir, "/results/pseudogene_ORs_distance_to_enhancers_Trim66_KO.pdf", sep=''), plot_step3.2_3_c, width = 50, height = 25, units ="cm")


#####################################################################################################################
#3.3 plot functional ORs and TAARs next to another using relative chromosome coordinates
###################################################################################################################### add one column to indicate if this receptor is significantly chaged in homo vs het
# add one column to show if the two ORs are above 1 Mb apart. Longzhi said :因为我画过or间距（log10之后）的histogram，1Mb那里是个谷
res_functional_OR_and_TAAR_step3.3 <- res_functional_OR_and_TAAR_step3.1
res_functional_OR_and_TAAR_step3.3$above_1Mb_between_consective_ORs <- diff(c(1, res_functional_OR_and_TAAR_step3.3$genomic_start)) >= 1000000 &
  diff(c(1, res_functional_OR_and_TAAR_step3.3$reversed_chr_number)) == 0 # add a row with number 1 for diff function. Also make sure it is in the same chromosome
# change the "above_1Mb_between_consective_ORs" value of the first and last of chr as TRUE
for (i in unique(res_functional_OR_and_TAAR_step3.3$chr_name)) {
  res_functional_OR_and_TAAR_step3.3[res_functional_OR_and_TAAR_step3.3$chr_name == i, ]$above_1Mb_between_consective_ORs[1] <- TRUE #the first as TRUE
  length_chr_temp_step3.3 <- length(res_functional_OR_and_TAAR_step3.3[res_functional_OR_and_TAAR_step3.3$chr_name == i, ]$above_1Mb_between_consective_ORs)
  res_functional_OR_and_TAAR_step3.3[res_functional_OR_and_TAAR_step3.3$chr_name == i, ]$above_1Mb_between_consective_ORs[length_chr_temp_step3.3] <- TRUE #the last as TRUE
}
# add another column as the x value for one-by-one plotting
res_functional_OR_and_TAAR_step3.3$chr_coordinate_for_plot <- NA
for (i in unique(res_functional_OR_and_TAAR_step3.3$chr_name)) {
  res_functional_OR_and_TAAR_step3.3[res_functional_OR_and_TAAR_step3.3$chr_name == i, ]$chr_coordinate_for_plot <- 1:sum(res_functional_OR_and_TAAR_step3.3$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}

res_functional_OR_and_TAAR_dataframe_step3.3 <- as.data.frame(res_functional_OR_and_TAAR_step3.3)
plot_step3.3 <- ggplot(data=res_functional_OR_and_TAAR_dataframe_step3.3, aes(x=chr_coordinate_for_plot, y=reversed_chr_number, 
                                                                              colour=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate,
                                                                              alpha=if_significant_upregulate_or_downregulate, shape=as.factor(genomic_strand))) +
  geom_point(size=3) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_shape_manual(values = c("\u25C4","\u25BA")) + #to choose other shapes like different triangles. https://stackoverflow.com/questions/30742379/creating-new-shape-palettes-in-ggplot2-and-other-r-graphics
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,1,1)) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("relative chr coordinate") + ylab("chr number") +
  scale_y_continuous(breaks = -1:-20, labels = c(1:19, "X"))  + # set the breaks to show as indicated, and also change the labels as indicated
  labs(colour = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "\u25BA plus strand\n\u25C4 minus strand") + # change legend title
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), #adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add columns for plotting the segments to separate OR clusters that are above 1Mb far
res_functional_OR_and_TAAR_dataframe_step3.3$segment_x <- ifelse(res_functional_OR_and_TAAR_dataframe_step3.3$above_1Mb_between_consective_ORs, 
                                                                 res_functional_OR_and_TAAR_dataframe_step3.3$chr_coordinate_for_plot-0.5, # subtract 0.5 so that the segments plotted will be on the left of triangle representing genes
                                                                 NA)
# change the "segemtn_x" of the last gene in chr to "chr_coordinate_for_plot"+0.5 other than "chr_coordinate_for_plot"-0.5
for (i in unique(res_functional_OR_and_TAAR_dataframe_step3.3$chr_name)) {
  res_functional_OR_and_TAAR_dataframe_step3.3[res_functional_OR_and_TAAR_dataframe_step3.3$chr_name == i, ]$segment_x[sum(res_functional_OR_and_TAAR_dataframe_step3.3$chr_name == i)] <- 
    res_functional_OR_and_TAAR_dataframe_step3.3[res_functional_OR_and_TAAR_dataframe_step3.3$chr_name == i, ]$chr_coordinate_for_plot[sum(res_functional_OR_and_TAAR_dataframe_step3.3$chr_name == i)] + 0.5
}
# add y_start and y_end for plot segment
res_functional_OR_and_TAAR_dataframe_step3.3$segment_y_start <- res_functional_OR_and_TAAR_dataframe_step3.3$reversed_chr_number-0.05
res_functional_OR_and_TAAR_dataframe_step3.3$segment_y_end <- res_functional_OR_and_TAAR_dataframe_step3.3$reversed_chr_number+0.05
# The geometric shape may look weird if you print to pdf file. check this website for solution:
# https://stackoverflow.com/questions/25053072/how-to-draw-half-filled-points-in-r-preferably-using-ggplot/34083063#34083063
grDevices::cairo_pdf(paste(resultDir, "/functional_ORs_TAARs_changed_in_genome_Trim66_KO.pdf", sep=""), width = 21, height = 27.9, family="Arial Unicode MS")
plot_step3.3 + annotate("segment", x = res_functional_OR_and_TAAR_dataframe_step3.3$segment_x, 
                        xend = res_functional_OR_and_TAAR_dataframe_step3.3$segment_x, 
                        y = res_functional_OR_and_TAAR_dataframe_step3.3$segment_y_start, 
                        yend = res_functional_OR_and_TAAR_dataframe_step3.3$segment_y_end,
                        size=0.3) # annotate to add the segment to separate the ORs if they are 1 Mb apart
dev.off()

#####################################################################################################################
#3.4 another way to plot pseudogene ORs next to another
###################################################################################################################### add one column to indicate if this receptor is significantly chaged in homo vs het
res_pseudogene_OR_step3.4 <- res_pseudogene_OR_step3.2
# add one column to show if the two ORs are above 1 Mb apart. Longzhi said :因为我画过or间距（log10之后）的histogram，1mb那里是个谷
res_pseudogene_OR_step3.4$above_1Mb_between_consective_ORs <- diff(c(1, res_pseudogene_OR_step3.4$genomic_start)) >= 1000000 &
  diff(c(1, res_pseudogene_OR_step3.4$reversed_chr_number)) == 0 # add a row with number 1 for diff function. Also make sure it is in the same chromosome
# change the "above_1Mb_between_consective_ORs" value of the first and last of chr as TRUE
for (i in unique(res_pseudogene_OR_step3.4$chr_name)) {
  res_pseudogene_OR_step3.4[res_pseudogene_OR_step3.4$chr_name == i, ]$above_1Mb_between_consective_ORs[1] <- TRUE #the first as TRUE
  length_chr_temp_step3.4 <- length(res_pseudogene_OR_step3.4[res_pseudogene_OR_step3.4$chr_name == i, ]$above_1Mb_between_consective_ORs)
  res_pseudogene_OR_step3.4[res_pseudogene_OR_step3.4$chr_name == i, ]$above_1Mb_between_consective_ORs[length_chr_temp_step3.4] <- TRUE #the last as TRUE
}
# add another column as the x value for one-by-one plotting
res_pseudogene_OR_step3.4$chr_coordinate_for_plot <- NA
for (i in unique(res_pseudogene_OR_step3.4$chr_name)) {
  res_pseudogene_OR_step3.4[res_pseudogene_OR_step3.4$chr_name == i, ]$chr_coordinate_for_plot <- 1:sum(res_pseudogene_OR_step3.4$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}

res_pseudogene_OR_dataframe_step3.4 <- as.data.frame(res_pseudogene_OR_step3.4)
plot_step3.4 <- ggplot(data=res_pseudogene_OR_dataframe_step3.4, aes(x=chr_coordinate_for_plot, y=reversed_chr_number, 
                                                                     colour=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate,
                                                                     alpha=if_significant_upregulate_or_downregulate, shape=as.factor(genomic_strand))) +
  geom_point(size=8) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_shape_manual(values = c("\u25C4","\u25BA")) + #to choose other shapes like different triangles. https://stackoverflow.com/questions/30742379/creating-new-shape-palettes-in-ggplot2-and-other-r-graphics
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,1,1)) +
  ggtitle('Pseudogenes ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("relative chr coordinate") + ylab("chr number") +
  scale_y_continuous(breaks = -1:-20, labels = c(1:19, "X"))  + # set the breaks to show as indicated, and also change the labels as indicated
  labs(colour = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "\u25BA plus strand\n\u25C4 minus strand") + # change legend title
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), #adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add columns for plotting the segments to separate OR clusters that are above 1Mb far
res_pseudogene_OR_dataframe_step3.4$segment_x <- ifelse(res_pseudogene_OR_dataframe_step3.4$above_1Mb_between_consective_ORs, 
                                                        res_pseudogene_OR_dataframe_step3.4$chr_coordinate_for_plot-0.5, 
                                                        NA)
# change the "segemtn_x" of the last gene in chr to "chr_coordinate_for_plot"+0.5 other than "chr_coordinate_for_plot"-0.5
for (i in unique(res_pseudogene_OR_dataframe_step3.4$chr_name)) {
  res_pseudogene_OR_dataframe_step3.4[res_pseudogene_OR_dataframe_step3.4$chr_name == i, ]$segment_x[sum(res_pseudogene_OR_dataframe_step3.4$chr_name == i)] <- 
    res_pseudogene_OR_dataframe_step3.4[res_pseudogene_OR_dataframe_step3.4$chr_name == i, ]$chr_coordinate_for_plot[sum(res_pseudogene_OR_dataframe_step3.4$chr_name == i)] + 0.5
}
# add y_start and y_end for plot segment
res_pseudogene_OR_dataframe_step3.4$segment_y_start <- res_pseudogene_OR_dataframe_step3.4$reversed_chr_number-0.15
res_pseudogene_OR_dataframe_step3.4$segment_y_end <- res_pseudogene_OR_dataframe_step3.4$reversed_chr_number+0.15
# The geometric shape may look weird if you print to pdf file. check this website for solution:
# https://stackoverflow.com/questions/25053072/how-to-draw-half-filled-points-in-r-preferably-using-ggplot/34083063#34083063
grDevices::cairo_pdf(paste(resultDir, "/pseudogenes_ORs_changed_in_genome_Trim66_KO.pdf", sep=""), width = 21, height = 27.9, family="Arial Unicode MS")
plot_step3.4 + annotate("segment", x = res_pseudogene_OR_dataframe_step3.4$segment_x, 
                        xend = res_pseudogene_OR_dataframe_step3.4$segment_x, 
                        y = res_pseudogene_OR_dataframe_step3.4$segment_y_start, 
                        yend = res_pseudogene_OR_dataframe_step3.4$segment_y_end,
                        size=0.8) # annotate to add the segment to separate the ORs if they are 1 Mb apart
dev.off()

#####################################################################################################################
#3.5 to check if significantly changed genes are in specific genome sites along real chr coordinates (plot all of the genes)
###################################################################################################################### add one column to indicate if this receptor is significantly chaged in homo vs het
res_all_step3.5 <- res_all
res_all_step3.5$if_significant_homoVShet <- res_all_step3.5$padj<0.05 & abs(res_all_step3.5$log2FoldChange)>=0.585
# replace NA with FALSE
res_all_step3.5[is.na(res_all_step3.5$padj), ]$if_significant_homoVShet <- FALSE
# add one column to indicate if this gene is significantly upregulated or downregulated in homo vs het
res_all_step3.5$if_significant_upregulate_or_downregulate <- "Not_changed"
# This is better. You can use for loop as the above code, but it takes longer as there are much more genes.
res_all_step3.5[which(res_all_step3.5$padj<0.05 & res_all_step3.5$log2FoldChange >= 0.585), ]$if_significant_upregulate_or_downregulate <- "Upregulated_significantly"
res_all_step3.5[which(res_all_step3.5$padj<0.05 & res_all_step3.5$log2FoldChange <= -0.585), ]$if_significant_upregulate_or_downregulate <- "Downregulated_significantly"
# add one column to inverse the chromosome numbers to negative numbers. Also remove chrY
res_all_remove_chrY_step3.5 <- res_all_step3.5[!grepl("^chrY", res_all_step3.5$chr_name), ]
res_all_remove_chrY_step3.5$reversed_chr_number <- ifelse(res_all_remove_chrY_step3.5$chr_name == "chrX", -20, as.numeric(substring(res_all_remove_chrY_step3.5$chr_name, 4)) * -1) #check substring function, https://stackoverflow.com/questions/17215789/extract-a-substring-according-to-a-pattern
# plot the ORs based on chromosome numbers, and strand direction
res_all_remove_chrY_dataframe_step3.5 <- as.data.frame(res_all_remove_chrY_step3.5)
plot_step3.5 <- ggplot(data=res_all_remove_chrY_dataframe_step3.5, aes(x=genomic_start, y=reversed_chr_number, 
                                                                       colour=if_significant_upregulate_or_downregulate, fill=if_significant_upregulate_or_downregulate,
                                                                       alpha=if_significant_upregulate_or_downregulate, shape=as.factor(genomic_strand))) +
  geom_point(size=2, position = position_jitter(width = 0, height = 0.4)) +  #"alpha" to set transparancy. The graphical argument used to specify point shapes is "pch", 16 is filled circle. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  scale_shape_manual(values = c("\u25C4","\u25BA")) + #to choose other shapes like different triangles. https://stackoverflow.com/questions/30742379/creating-new-shape-palettes-in-ggplot2-and-other-r-graphics
  scale_color_manual(values = c("springgreen3", "gray90", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray90", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25))) +
  ggtitle('All genes including pseudogenes changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("chr coordinate") + ylab("chr number") +
  xlim(0, 2e8) + # 2.0e+8 as the maximum of plot_step3.5
  scale_y_continuous(breaks = -1:-20, labels = c(1:19, "X"))  + # set the breaks to show as indicated, and also change the labels as indicated
  labs(colour = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "\u25BA plus strand\n\u25C4 minus strand") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
# The geometric shape may look weird if you print to pdf file. check this website for solution:
# https://stackoverflow.com/questions/25053072/how-to-draw-half-filled-points-in-r-preferably-using-ggplot/34083063#34083063
grDevices::cairo_pdf(paste(resultDir, "/all_genes_changed_in_genome_Trim66_KO_whole_genome.pdf", sep=""), width = 21, height = 27.9, family="Arial Unicode MS")
plot_step3.5
dev.off()
# compare plot_step3.5 with plot_step3.1 and plot_step3.2
arrange_plots_step3.5 <- ggarrange(plot_step3.5, plot_step3.1, plot_step3.2, ncol = 3, nrow = 1)
ggsave(paste(baseDir, "/results/combined_genes_changed_in_genome_Trim66_KO.pdf", sep=''), arrange_plots_step3.5, width = 120, height = 60, units ="cm")


#####################################################################################################################
          # Plot Step 4: plot OR and TAAR changes with relative chromosome information
#####################################################################################################################
# genes in res_all are sorted by chr and chr coordinates, which is good for the following analysis
# 4.1. plot functional OR and TAAR genes along relative chromosomes
# add one column as linear x axis for plotting
res_functional_OR_and_TAAR_dataframe_step4.1 <- res_functional_OR_and_TAAR_dataframe_step3.1
res_functional_OR_and_TAAR_dataframe_step4.1$x_axis_for_plot_relative_chr_position <- 1:nrow(res_functional_OR_and_TAAR_dataframe_step4.1)
# add one column to indicate TAAR or OR
res_functional_OR_and_TAAR_dataframe_step4.1$TAAR_or_OR <- ifelse(grepl("^Taar", rownames(res_functional_OR_and_TAAR_dataframe_step4.1)), "TAAR", "OR")
# add one column to change NA in log2FoldChange to 0 so that those data will be not removed when plotting
res_functional_OR_and_TAAR_dataframe_step4.1$log2FoldChange_NA_to_0 <- res_functional_OR_and_TAAR_dataframe_step4.1$log2FoldChange
res_functional_OR_and_TAAR_dataframe_step4.1[is.na(res_functional_OR_and_TAAR_dataframe_step4.1$log2FoldChange),]$log2FoldChange_NA_to_0 <- 0
#plot
# generate a universal variable for chr colors (all 21 chr including chrY. chrX and chrY use different colors) used for all plot.
## not using scale_color_brewer as there are 21 different chr, but most color palettes have 8-12 colors. check http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html  and https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21) 
col_for_chr_step4.1 <- col_for_chr_universal
plot_step4.1 <- ggplot(data=res_functional_OR_and_TAAR_dataframe_step4.1, aes(x=x_axis_for_plot_relative_chr_position, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR), #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 alpha=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR)), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21, stroke=0) + #stroke to control point border thickness
  geom_line(aes(x=x_axis_for_plot_relative_chr_position, y=8,color=factor(chr_name, levels = gtools::mixedsort(unique(chr_name)))), #make chr_name as factor and relevel according the order chr1, chr2, chr3.... instead of the default order of chr1, chr10, chr11... check https://www.biostars.org/p/9504853/#9504858
            linetype="solid", size=20) + #plot the chr information with different chr using different colors
  scale_color_manual(values = col_for_chr_step4.1,
                     guide = guide_legend(nrow=2)) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_alpha_manual(values = c(1,0.5,1,1,0.8,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("segment", x = res_functional_OR_and_TAAR_dataframe_step4.1[res_functional_OR_and_TAAR_dataframe_step4.1$chr_name == "chr5", ]$x_axis_for_plot_relative_chr_position, 
           xend = res_functional_OR_and_TAAR_dataframe_step4.1[res_functional_OR_and_TAAR_dataframe_step4.1$chr_name == "chr5", ]$x_axis_for_plot_relative_chr_position, 
           y = 7.73, 
           yend = 8.27,
           size=1, # annotate to add the segment as only 1 OR in chr5 so the geom_line cannot draw line for chr5
           color= col_for_chr_step4.1[5]) + # color is the 5th of the color palette generated above
  ylim(-8, 8.3) + # range(res_functional_OR_and_TAAR_dataframe_step4.1$log2FoldChange_NA_to_0) to check the y axis range
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Relative chromosome coordinates") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "chr name: note that \nno functional OR genes \nin chr12, chr18, and chrY. \nThis is just legend. \ndelete them in legend when making real figures", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

ggsave(paste(baseDir, "/results/functional_ORs_TAARs_relative_chr_Trim66_KO.pdf", sep=''), plot_step4.1, width = 120, height = 60, units ="cm")

# 4.2. plot pseudogene OR genes along relative chromosomes
# add one column as linear x axis for plotting
res_pseudogene_OR_dataframe_step4.2 <- res_pseudogene_OR_dataframe_step3.2
res_pseudogene_OR_dataframe_step4.2$x_axis_for_plot_relative_chr_position <- 1:nrow(res_pseudogene_OR_dataframe_step4.2)
# add one column to change NA in log2FoldChange to 0 so that those data will be not removed when plotting
res_pseudogene_OR_dataframe_step4.2$log2FoldChange_NA_to_0 <- res_pseudogene_OR_dataframe_step4.2$log2FoldChange
res_pseudogene_OR_dataframe_step4.2[is.na(res_pseudogene_OR_dataframe_step4.2$log2FoldChange),]$log2FoldChange_NA_to_0 <- 0
#plot
col_for_chr_step4.2 <- col_for_chr_universal
plot_step4.2 <- ggplot(data=res_pseudogene_OR_dataframe_step4.2, aes(x=x_axis_for_plot_relative_chr_position, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21, stroke=0) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  geom_line(aes(x=x_axis_for_plot_relative_chr_position, y=7,color=factor(chr_name, levels = gtools::mixedsort(unique(chr_name)))), #make chr_name as factor and relevel according the order chr1, chr2, chr3.... instead of the default order of chr1, chr10, chr11...
            linetype="solid", size=20) + #plot the chr information with different chr using different colors
  scale_color_manual(values = col_for_chr_step4.2,
                     guide = guide_legend(nrow=2)) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("segment", x = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chr12", ]$x_axis_for_plot_relative_chr_position, 
           xend = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chr12", ]$x_axis_for_plot_relative_chr_position, 
           y = 6.78, 
           yend = 7.22,
           size=3, # annotate to add the segment as only 1 pseudogene OR in chr12,13,15,and chrX so the geom_line cannot draw line for those chr. Note that there are no pseugogenes in chr8, chr18, chrY.
           color=col_for_chr_step4.2[12]) + # for chr12
  annotate("segment", x = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chr13", ]$x_axis_for_plot_relative_chr_position, 
           xend = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chr13", ]$x_axis_for_plot_relative_chr_position, 
           y = 6.78, 
           yend = 7.22,
           size=3, # annotate to add the segment as only 1 pseudogene OR in chr12,13,15,and chrX so the geom_line cannot draw line for them
           color=col_for_chr_step4.2[13]) + # for chr13
  annotate("segment", x = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chr15", ]$x_axis_for_plot_relative_chr_position, 
           xend = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chr15", ]$x_axis_for_plot_relative_chr_position, 
           y = 6.78, 
           yend = 7.22,
           size=3, # annotate to add the segment as only 1 pseudogene OR in chr12,13,15,and chrX so the geom_line cannot draw line for them
           color=col_for_chr_step4.2[15]) + # for chr15
  annotate("segment", x = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chrX", ]$x_axis_for_plot_relative_chr_position, 
           xend = res_pseudogene_OR_dataframe_step4.2[res_pseudogene_OR_dataframe_step4.2$chr_name == "chrX", ]$x_axis_for_plot_relative_chr_position, 
           y = 6.78, 
           yend = 7.22,
           size=3, # annotate to add the segment as only 1 pseudogene OR in chr12,13,15,and chrX so the geom_line cannot draw line for them
           color=col_for_chr_step4.2[20]) + # chrX is the 20th
  ylim(-6, 7.5) + # range(res_pseudogene_OR_dataframe_step4.2$log2FoldChange_NA_to_0) to check the y axis range
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Relative chromosome coordinates") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "chr name: note that \nno pseudogene OR genes \nin chr8, chr18, and chrY. \nThis is just legend. \ndelete them in legend when making real figures", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

ggsave(paste(baseDir, "/results/pseudogene_ORs_relative_chr_Trim66_KO.pdf", sep=''), plot_step4.2, width = 120, height = 60, units ="cm")

# Step 5: TPM calculation -------------------------------------------------

#####################################################################################################################
      # Plot Step 5: calculate TPM value from counts
#####################################################################################################################
#    this code is used to generate  TPM values from featureCounts data
#     https://www.jianshu.com/p/eece90bdddf9
# countdata with the gene length are obtained from featureCounts data, see RNA_seq_code.sh step 6 for details.
countdata_step5 <- read.delim(file.path(dataDir, "all_trimmomatic_featureCounts_for_TPM_calculation_with_header.txt"), header=T, sep="\t",
                              row.names=1, as.is=T)
countdata_good_genes_and_pseudogenes_step5 <- countdata_step5[match(original_counts$gene_id_Ensembl, rownames(countdata_step5)), ]
match(rownames(countdata_good_genes_and_pseudogenes_step5), original_counts$gene_id_Ensembl)
rownames(countdata_good_genes_and_pseudogenes_step5) <- rownames(original_counts)
original_colnames_step5 <- colnames(countdata_good_genes_and_pseudogenes_step5)
colnames(countdata_good_genes_and_pseudogenes_step5) <- c(original_colnames_step5[1], paste(original_colnames_step5[-1], "_raw_counts", sep = ""))
prefix_step5<-"Trim66_KO_het_vs_homo_RNA-seq_TPM"#设置输出文件前缀名
#-----TPM Calculation------
kb_step5 <- countdata_good_genes_and_pseudogenes_step5$gene_length / 1000
rpk_step5 <- countdata_good_genes_and_pseudogenes_step5[,-1] / kb_step5 # the first column of countdata_good_genes_and_pseudogenes_step5 is gene_length
tpm_step5 <- t(t(rpk_step5)/colSums(rpk_step5) * 1000000)
colnames(tpm_step5) <- paste(original_colnames_step5[-1], "_TPM", sep = "")
TPM_final_for_save_step5 <- cbind(tpm_step5, original_counts[,11:ncol(original_counts)], countdata_good_genes_and_pseudogenes_step5)
write.csv(TPM_final_for_save_step5,file=file.path(resultDir, paste0(prefix_step5,".csv")))


# Step 6: plot receptor changes against TPM -------------------------------

#####################################################################################################################
      # Plot Step 6: plot OR and TAAR changes with their expression level to see if they are correlated
#####################################################################################################################
#6.1 functional ORs and TAARs
all(rownames(TPM_final_for_save_step5) == rownames(res_all)) #res_all and TPM_final_for_save_step5 have same rownames in order
TPM_functional_OR_and_TAAR_step6.1 <- TPM_final_for_save_step5[grepl("Olfr|Taar[2-9]", rownames(res_all)) & grepl("protein_coding", res_all$gene_type), ]
all(rownames(TPM_functional_OR_and_TAAR_step6.1) == rownames(res_functional_OR_and_TAAR_dataframe_step6.1))
# plot the receptor changes against expression levels (average TPM from het)
res_functional_OR_and_TAAR_dataframe_step6.1 <- res_functional_OR_and_TAAR_dataframe_step4.1
res_functional_OR_and_TAAR_dataframe_step6.1$TPM_average_het_5samples <- rowMeans(TPM_functional_OR_and_TAAR_step6.1[, 1:5])
plot_step6.1_1 <- ggplot(data=res_functional_OR_and_TAAR_dataframe_step6.1, aes(x=TPM_average_het_5samples, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR), #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR),
                 alpha=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR)), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_alpha_manual(values = c(1,0.5,1,1,0.8,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 100, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_functional_OR_and_TAAR_dataframe_step6.1$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step6.1$TPM_average_het_5samples), 4),
                         ", p value: ",
                         round(cor.test(res_functional_OR_and_TAAR_dataframe_step6.1$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step6.1$TPM_average_het_5samples)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("TPM of olfactory receptor genes") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
res_functional_OR_and_TAAR_dataframe_step6.1$name <- rownames(res_functional_OR_and_TAAR_dataframe_step6.1)
plot_step6.1_2 <- ggplot(data=res_functional_OR_and_TAAR_dataframe_step6.1, aes(x=TPM_average_het_5samples, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR), #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR),
                 alpha=interaction(if_significant_upregulate_or_downregulate, TAAR_or_OR)), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1", "blue1", "black", "magenta1")) +
  scale_alpha_manual(values = c(1,0.5,1,1,0.8,1)) +
  geom_text(data=subset(res_functional_OR_and_TAAR_dataframe_step6.1, (padj<0.05 & abs(log2FoldChange)>=0.585) ==TRUE), 
            aes(label = name), 
            size =5, colour='blue') +  # label=name in the aes() suggest the text to label
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 100, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_functional_OR_and_TAAR_dataframe_step6.1$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step6.1$TPM_average_het_5samples), 4),
                         ", p value: ",
                         round(cor.test(res_functional_OR_and_TAAR_dataframe_step6.1$log2FoldChange_NA_to_0, res_functional_OR_and_TAAR_dataframe_step6.1$TPM_average_het_5samples)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional ORs and TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("TPM of olfactory receptor genes") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

arrange_step6.1 <- ggarrange(plot_step6.1_1, plot_step6.1_2, ncol = 1, nrow = 2)
ggsave(paste(baseDir, "/results/functional_ORs_TAARs_changes_against_TPM_Trim66_KO.pdf", sep=''), arrange_step6.1, width = 120, height = 60, units ="cm")

#6.2 pseudogene ORs
TPM_pseudogene_OR_step6.2 <- TPM_final_for_save_step5[grepl('^Olfr', rownames(res_all)) & grepl('pseudogene', res_all$gene_type), ]
res_pseudogene_OR_dataframe_step6.2 <- res_pseudogene_OR_dataframe_step4.2
all(rownames(TPM_pseudogene_OR_step6.2) == rownames(res_pseudogene_OR_dataframe_step6.2))
# plot the receptor changes against expression levels (average TPM from het)
res_pseudogene_OR_dataframe_step6.2$TPM_average_het_5samples <- rowMeans(TPM_pseudogene_OR_step6.2[, 1:5])
plot_step6.2_1 <- ggplot(data=res_pseudogene_OR_dataframe_step6.2, aes(x=TPM_average_het_5samples, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1,1,0.8,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 20, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_pseudogene_OR_dataframe_step6.2$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step6.2$TPM_average_het_5samples), 4),
                         ", p value: ",
                         round(cor.test(res_pseudogene_OR_dataframe_step6.2$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step6.2$TPM_average_het_5samples)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("TPM of olfactory receptor genes") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
res_pseudogene_OR_dataframe_step6.2$name <- rownames(res_pseudogene_OR_dataframe_step6.2)
plot_step6.2_2 <- ggplot(data=res_pseudogene_OR_dataframe_step6.2, aes(x=TPM_average_het_5samples, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1,1,0.8,1)) +
  geom_text(data=subset(res_pseudogene_OR_dataframe_step6.2, (padj<0.05 & abs(log2FoldChange)>=0.585) ==TRUE), 
            aes(label = name), 
            size =5, colour='blue') +  # label=name in the aes() suggest the text to label
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 20, y = 6, 
           label = paste("Pearson Correlation: ", 
                         round(cor(res_pseudogene_OR_dataframe_step6.2$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step6.2$TPM_average_het_5samples), 4),
                         ", p value: ",
                         round(cor.test(res_pseudogene_OR_dataframe_step6.2$log2FoldChange_NA_to_0, res_pseudogene_OR_dataframe_step6.2$TPM_average_het_5samples)$p.value, 4), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("TPM of olfactory receptor genes") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant", shape = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

arrange_plots_step6.2 <- ggarrange(plot_step6.2_1, plot_step6.2_2, ncol = 1, nrow = 2)
ggsave(paste(baseDir, "/results/pseudogene_ORs_changes_against_TPM_Trim66_KO.pdf", sep=''), arrange_plots_step6.2, width = 120, height = 60, units ="cm")


# Step 7: plot changes of cell markers, TFs and others --------------------

#####################################################################################################################
  # Plot Step 7: plot cell markers and other interested genes
#####################################################################################################################
# 7.1 (1) marker genes with normalized TPM
HBC_marker_genes_step7.1 <-c("Krt5","Krt14","Trp63") 
GBC_marker_genes_step7.1 <- c("Hes6","Sox2","Kit","Ascl1","Neurog1","Neurod1")
immature_OSN_marker_genes_step7.1 <- c("Gap43","Gnas","Gng8","Hdac2","Stmn4")
mature_OSN_marker_genes_step7.1 <- c("Omp","Gnal","Gng13","Adcy3","Cnga4","Cnga2","Cngb1","Slc17a6","Stoml3")
Ms4a_cell_marker_genes_step7.1 <- c("Car2","Cnga3")
microvillous_cell_marker_genes_step7.1 <- c("Trpm5", "Pou2f3") #from paper: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07528-y
sustentacular_cell_marker_genes_step7.1 <- c("Cyp2g1","Cyp1a2","Hes1")
bowman_gland_marker_genes_step7.1 <- c("Sox9", "Sox10")
# combine all the marker genes
all_marker_genes_step7.1 <- c(HBC_marker_genes_step7.1, GBC_marker_genes_step7.1, immature_OSN_marker_genes_step7.1, 
                              mature_OSN_marker_genes_step7.1, Ms4a_cell_marker_genes_step7.1, microvillous_cell_marker_genes_step7.1,
                              sustentacular_cell_marker_genes_step7.1, bowman_gland_marker_genes_step7.1)
marker_genes_TPM_data_step7.1 <- TPM_final_for_save_step5[match(all_marker_genes_step7.1, rownames(TPM_final_for_save_step5)), 1:10]
marker_genes_TPM_normalize_het_to_1_step7.1 <- marker_genes_TPM_data_step7.1 / rowMeans(marker_genes_TPM_data_step7.1[,1:5]) #normalize to the mean of het samples
marker_genes_TPM_normalize_het_to_1_step7.1$name <- rownames(marker_genes_TPM_normalize_het_to_1_step7.1)
marker_genes_TPM_normalize_het_to_1_step7.1 <- melt(marker_genes_TPM_normalize_het_to_1_step7.1, id.vars = "name") #melt function from reshape2 package to convert dataframe to the format easier for ggplot
marker_genes_TPM_normalize_het_to_1_step7.1$genotype <- "het"
marker_genes_TPM_normalize_het_to_1_step7.1[grep("homo", marker_genes_TPM_normalize_het_to_1_step7.1$variable), ]$genotype <- "homo"
marker_genes_TPM_normalize_het_to_1_step7.1$genotype <- factor(marker_genes_TPM_normalize_het_to_1_step7.1$genotype, levels=c("het", "homo")) # force ggplot to plot het first
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
marker_genes_TPM_normalize_het_to_1_step7.1$name <- factor(marker_genes_TPM_normalize_het_to_1_step7.1$name, levels = all_marker_genes_step7.1)
# follow this link for plotting median+quartiles, grouped plots: http://www.sthda.com/english/articles/32-r-graphics-essentials/132-plot-grouped-data-box-plot-bar-plot-and-more/
# if you want to plot mean+sem, follow this website: https://stackoverflow.com/questions/28819135/combined-bar-plot-and-points-in-ggplot2
plot_step7.1_1 <- ggpubr::ggbarplot(marker_genes_TPM_normalize_het_to_1_step7.1, x = "name", y = "value", add = c("mean_se"),
                                    fill = "gray95", size = 0.25, # size to change the size of outlines. 
                                    add.params = list(size=.25), # Add.params to change the size for the argument "add", the error bar in this case
                                    title = "Marker gene expression (TPM normalized to het) in Trim66 KO by ggpubr::ggbarplot", x.text.angle = 60, # rotate the x label
                                    ylab = "Fold change of TPM (mean+-sem)", xlab = "Marker genes",
                                    color = "genotype", palette = c("springgreen3","red1"),
                                    position = position_dodge(0.85)) + # This is better. Use it to avoid add("mean_se", "jitter") where you can not adjust points and error bar together, see the example in TAAR silencer analysis code. check https://github.com/kassambara/ggpubr/issues/130
  geom_point(aes(fill=genotype),
             size=2, shape=21, stroke = 0.2, # stroke control border thickness
             color="black",
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8)) +
  scale_fill_manual(values = rep("white", 2))
# obtain the DESeq2 results of marker genes
marker_genes_statistics_data_step7.1 <- res_all[match(all_marker_genes_step7.1, rownames(res_all)), c(1:2,5:6)]
marker_genes_statistics_data_step7.1$TPM_mean_het_5samples <- rowMeans(marker_genes_TPM_data_step7.1[,1:5])
marker_genes_statistics_data_step7.1$TPM_mean_homo_5samples <- rowMeans(marker_genes_TPM_data_step7.1[,6:10])
# Draw ggplot2 plot with text only with p value information
plot_step7.1_2 <- ggplot() +  #check https://statisticsglobe.com/plot-only-text-in-r
  annotate("text",
           x = 0,
           y = 0,
           size = 3, #adjust the size if the data are not fully printed
           label = paste(capture.output(as.data.frame(marker_genes_statistics_data_step7.1)), collapse = "\n")) + # check https://stackoverflow.com/questions/14326573/print-a-data-frame-in-the-white-space-of-a-plot
  theme_void()

#####################################################################################
          #note: using ggplot2 is better
   #Next time use this code to plot mean+-se with points. Get the same plot as using ggpubr::ggbarplot
#####################################################################################
# 7.1 (2)marker genes using normalized TPM with ggplot
#### This is better. Next time use ggplot2 to make this kind of figures.
marker_genes_TPM_normalize_means_se_step7.1_2 <- marker_genes_TPM_normalize_het_to_1_step7.1 %>%  group_by(name, genotype) %>% 
  summarize(mean_TPM=mean(value), 
            sd_TPM=sd(value), 
            Numbers=n(), 
            se=sd_TPM/sqrt(Numbers), 
            upper_limit=mean_TPM+se, 
            lower_limit=mean_TPM-se) 
marker_genes_TPM_normalize_means_se_step7.1_2
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
marker_genes_TPM_normalize_means_se_step7.1_2$name <- factor(marker_genes_TPM_normalize_means_se_step7.1_2$name, levels = all_marker_genes_step7.1)
plot_step7.1_3 <- ggplot() + 
  geom_bar(marker_genes_TPM_normalize_means_se_step7.1_2,
           mapping=aes(x=name, y=mean_TPM, color=genotype),
           stat="identity", 
           position = position_dodge(0.8),
           width = 0.7, #adjust the bar width
           size=0.25, # adjust the line thickness
           fill="gray95") + 
  geom_errorbar(marker_genes_TPM_normalize_means_se_step7.1_2,
                mapping=aes(x=name, ymin=lower_limit, ymax=upper_limit, color=genotype),
                position = position_dodge(0.8), width = .2, size=0.25) +
  scale_color_manual(values = c("springgreen3", "red1")) +
  geom_point(marker_genes_TPM_normalize_het_to_1_step7.1,
             mapping=aes(x=name, y=value, group=genotype, fill=genotype),
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8),
             size=2, shape=21, stroke = 0.2, # stroke control border thickness
             color="black") +
  scale_fill_manual(values = rep("white", 2)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('Marker gene expression (TPM normalized to het) in Trim66 KO by ggplot') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Marker genes") + ylab("Fold change of TPM (mean+-sem)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

arrange_plots_step7.1 <- ggarrange(plot_step7.1_1, plot_step7.1_3, plot_step7.1_2, ncol = 1, nrow = 3)
ggsave(paste(baseDir, "/results/marker_genes_normalized_TPM_Trim66_KO.pdf", sep=''), arrange_plots_step7.1, width = 50, height = 50, units ="cm")

# 7.2 (1)important transcription factors normalized TPM
key_TF_genes_step7.2 <- c("Lhx2","Emx2","Ebf1", "Ebf2", "Ebf3", "Ebf4", "Atf5", "Kdm1a", "Ldb1", "Bptf") #from Stavros' Nature, Cell, and eLife paper on enhancer and Lsd1
key_TF_genes_TPM_data_step7.2_1 <- TPM_final_for_save_step5[match(key_TF_genes_step7.2, rownames(TPM_final_for_save_step5)), 1:10]
key_TF_genes_TPM_normalize_het_to_1_step7.2_1 <- key_TF_genes_TPM_data_step7.2_1 / rowMeans(key_TF_genes_TPM_data_step7.2_1[,1:5]) #normalize to the mean of het samples
key_TF_genes_TPM_normalize_het_to_1_step7.2_1$name <- rownames(key_TF_genes_TPM_normalize_het_to_1_step7.2_1)
key_TF_genes_TPM_normalize_het_to_1_step7.2_1 <- melt(key_TF_genes_TPM_normalize_het_to_1_step7.2_1, id.vars = "name") #melt function from reshape2 package to convert dataframe to the format easier for ggplot
key_TF_genes_TPM_normalize_het_to_1_step7.2_1$genotype <- "het"
key_TF_genes_TPM_normalize_het_to_1_step7.2_1[grep("homo", key_TF_genes_TPM_normalize_het_to_1_step7.2_1$variable), ]$genotype <- "homo"
key_TF_genes_TPM_normalize_het_to_1_step7.2_1$genotype <- factor(key_TF_genes_TPM_normalize_het_to_1_step7.2_1$genotype, levels=c("het", "homo")) # force ggplot to plot het first
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
key_TF_genes_TPM_normalize_het_to_1_step7.2_1$name <- factor(key_TF_genes_TPM_normalize_het_to_1_step7.2_1$name, levels = key_TF_genes_step7.2)
# follow this link for plotting median+quartiles, grouped plots: http://www.sthda.com/english/articles/32-r-graphics-essentials/132-plot-grouped-data-box-plot-bar-plot-and-more/
# if you want to plot mean+sem, follow this website: https://stackoverflow.com/questions/28819135/combined-bar-plot-and-points-in-ggplot2
plot_step7.2_1_1 <- ggpubr::ggbarplot(key_TF_genes_TPM_normalize_het_to_1_step7.2_1, x = "name", y = "value", add = c("mean_se"),
                                      fill = "gray95", size = 0.25, # size to change the size of outlines. 
                                      add.params = list(size=.25), # Add.params to change the size for the argument "add", the error bar in this case
                                      title = "Transcription regulatory gene expression (TPM normalized to het) in Trim66 KO by ggpubr::ggbarplot", x.text.angle = 60, # rotate the x label
                                      ylab = "Fold change of TPM (mean+-sem)", xlab = "Transcription regulatory genes",
                                      color = "genotype", palette = c("springgreen3", "red1"),
                                      position = position_dodge(0.85)) + # This is better. Use it to avoid add("mean_se", "jitter") where you can not adjust points and error bar together, see the example in TAAR silencer analysis code. check https://github.com/kassambara/ggpubr/issues/130
  geom_point(aes(fill=genotype),
             size=4, shape=21, stroke = 0.2, # stroke control border thickness
             color="black",
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8)) +
  scale_fill_manual(values = rep("white", 2))
# obtain the DESeq2 results of TF genes
key_TF_genes_statistics_data_step7.2_1 <- res_all[match(key_TF_genes_step7.2, rownames(res_all)), c(1:2,5:6)]
key_TF_genes_statistics_data_step7.2_1$TPM_mean_het_5samples <- rowMeans(key_TF_genes_TPM_data_step7.2_1[,1:5])
key_TF_genes_statistics_data_step7.2_1$TPM_mean_homo_5samples <- rowMeans(key_TF_genes_TPM_data_step7.2_1[,6:10])
# Draw ggplot2 plot with text only with p value information
plot_step7.2_1_2 <- ggplot() +  #check https://statisticsglobe.com/plot-only-text-in-r
  annotate("text",
           x = 0,
           y = 0,
           size = 3, #adjust the size if the data are not fully printed
           label = paste(capture.output(as.data.frame(key_TF_genes_statistics_data_step7.2_1)), collapse = "\n")) + # check https://stackoverflow.com/questions/14326573/print-a-data-frame-in-the-white-space-of-a-plot
  theme_void()

#####################################################################################
#note: using ggplot2 is better
#####################################################################################
# 7.2 (2)marker genes using normalized TPM with ggplot
#### This is better. Next time use ggplot2 to make this kind of figures.
key_TF_genes_TPM_normalize_means_se_step7.2_2 <- key_TF_genes_TPM_normalize_het_to_1_step7.2_1 %>%  group_by(name, genotype) %>% 
  summarize(mean_TPM=mean(value), 
            sd_TPM=sd(value), 
            Numbers=n(), 
            se=sd_TPM/sqrt(Numbers), 
            upper_limit=mean_TPM+se, 
            lower_limit=mean_TPM-se) 
key_TF_genes_TPM_normalize_means_se_step7.2_2
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
key_TF_genes_TPM_normalize_means_se_step7.2_2$name <- factor(key_TF_genes_TPM_normalize_means_se_step7.2_2$name, levels = key_TF_genes_step7.2)
plot_step7.2_2 <- ggplot() + 
  geom_bar(key_TF_genes_TPM_normalize_means_se_step7.2_2,
           mapping=aes(x=name, y=mean_TPM, color=genotype),
           stat="identity", 
           position = position_dodge(0.8),
           width = 0.7, #adjust the bar width
           size=0.25, # adjust the line thickness
           fill="gray95") + 
  geom_errorbar(key_TF_genes_TPM_normalize_means_se_step7.2_2,
                mapping=aes(x=name, ymin=lower_limit, ymax=upper_limit, color=genotype),
                position = position_dodge(0.8), width = .2, size=0.25) +
  scale_color_manual(values = c("springgreen3", "red1")) +
  geom_point(key_TF_genes_TPM_normalize_het_to_1_step7.2_1,
             mapping=aes(x=name, y=value, group=genotype, fill=genotype),
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8),
             size=4, shape=21, stroke = 0.2, # stroke control border thickness
             color="black") +
  scale_fill_manual(values = rep("white", 2)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('Transcription regulatory gene expression (TPM normalized to het) in Trim66 KO by ggplot') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Transcription regulatory genes") + ylab("Fold change of TPM (mean+-sem)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

arrange_plots_step7.2_1and2 <- ggarrange(plot_step7.2_1_1, plot_step7.2_2, plot_step7.2_1_2, ncol = 1, nrow = 3)
ggsave(paste(baseDir, "/results/key_TF_genes_normalized_TPM_Trim66_KO.pdf", sep=''), arrange_plots_step7.2_1and2, width = 50, height = 50, units ="cm")


# 7.2 (3)important transcriptional factors original TPM
key_TF_genes_step7.2 <- c("Lhx2","Emx2","Ebf1", "Ebf2", "Ebf3", "Ebf4", "Atf5", "Kdm1a", "Ldb1", "Bptf") #from Stavros' Nature, Cell, and eLife paper on enhancer and Lsd1
key_TF_genes_TPM_data_step7.2_3 <- TPM_final_for_save_step5[match(key_TF_genes_step7.2, rownames(TPM_final_for_save_step5)), 1:10]
# statstical data of gene changes
key_TF_genes_statistics_data_step7.2_3 <- res_all[match(key_TF_genes_step7.2, rownames(res_all)), c(1:2,5:6)]
key_TF_genes_statistics_data_step7.2_3$TPM_mean_het_5samples <- rowMeans(key_TF_genes_TPM_data_step7.2_3[,1:5])
key_TF_genes_statistics_data_step7.2_3$TPM_mean_homo_5samples <- rowMeans(key_TF_genes_TPM_data_step7.2_3[,6:10])
# arrange data for plotting
key_TF_genes_TPM_data_step7.2_3$name <- rownames(key_TF_genes_TPM_data_step7.2_3)
key_TF_genes_TPM_data_step7.2_3 <- melt(key_TF_genes_TPM_data_step7.2_3, id.vars = "name") #melt function from reshape2 package to convert dataframe to the format easier for ggplot
key_TF_genes_TPM_data_step7.2_3$genotype <- "het"
key_TF_genes_TPM_data_step7.2_3[grep("homo", key_TF_genes_TPM_data_step7.2_3$variable), ]$genotype <- "homo"
key_TF_genes_TPM_data_step7.2_3$genotype <- factor(key_TF_genes_TPM_data_step7.2_3$genotype, levels=c("het", "homo")) # force ggplot to plot het first
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
key_TF_genes_TPM_data_step7.2_3$name <- factor(key_TF_genes_TPM_data_step7.2_3$name, levels = key_TF_genes_step7.2)
#### use a different method to add text here, as using the above step 7.2 (1) method to print text in plot and ggarrange does not work for broken y axis with ggbreak::scale_y_break function (the broken y axis will be back to normal by ggarrange)
plot_step7.2_3 <- ggpubr::ggbarplot(key_TF_genes_TPM_data_step7.2_3, x = "name", y = "value", add = c("mean_se"),
                                    fill = "gray95", size = 0.25, # size to change the size of outlines. 
                                    add.params = list(size=.25), # Add.params to change the size for the argument "add", the error bar and point in this case
                                    title = "Transcription regulatory gene expression (original TPM) \nin Trim66 KO by ggpubr::ggbarplot", x.text.angle = 60, # rotate the x label
                                    ylab = "TPM (mean+-sem)", xlab = "Transcription regulatory genes",
                                    color = "genotype", palette = c("springgreen3", "red1"),
                                    position = position_dodge(0.85)) + # This is better. Use it to avoid add("mean_se", "jitter"), see the example in TAAR silencer analysis code. check https://github.com/kassambara/ggpubr/issues/130
  geom_point(aes(fill=genotype),
             size=2, shape=21, stroke = 0.2, # stroke control border thickness
             color="black",
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8)) +
  scale_fill_manual(values = rep("white", 2)) +
  ggbreak::scale_y_break(c(300, 2000)) +# check https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html#introduction
  annotate("text", x = c("Ebf1"), y = 2400, 
           label = paste(capture.output(as.data.frame(key_TF_genes_statistics_data_step7.2_3)), collapse = "\n"), # check https://stackoverflow.com/questions/14326573/print-a-data-frame-in-the-white-space-of-a-plot
           size=1.5, colour = "black") # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console

#####################################################################################
#note: using ggplot2 is better
#####################################################################################
# 7.2 (4)marker genes using original raw TPM with ggplot
#### This is better. Next time use ggplot2 to make this kind of figures.
key_TF_genes_TPM_means_se_step7.2_4 <- key_TF_genes_TPM_data_step7.2_3 %>%  group_by(name, genotype) %>% 
  summarize(mean_TPM=mean(value), 
            sd_TPM=sd(value), 
            Numbers=n(), 
            se=sd_TPM/sqrt(Numbers), 
            upper_limit=mean_TPM+se, 
            lower_limit=mean_TPM-se) 
key_TF_genes_TPM_means_se_step7.2_4
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
key_TF_genes_TPM_means_se_step7.2_4$name <- factor(key_TF_genes_TPM_means_se_step7.2_4$name, levels = key_TF_genes_step7.2)
plot_step7.2_4 <- ggplot() + 
  geom_bar(key_TF_genes_TPM_means_se_step7.2_4,
           mapping=aes(x=name, y=mean_TPM, color=genotype),
           stat="identity", 
           position = position_dodge(0.8),
           width = 0.7, #adjust the bar width
           size=0.25, # adjust the line thickness
           fill="gray95") + 
  geom_errorbar(key_TF_genes_TPM_means_se_step7.2_4,
                mapping=aes(x=name, ymin=lower_limit, ymax=upper_limit, color=genotype),
                position = position_dodge(0.8), width = .2, size=0.25) +
  scale_color_manual(values = c("springgreen3", "red1")) +
  geom_point(key_TF_genes_TPM_data_step7.2_3,
             mapping=aes(x=name, y=value, group=genotype, fill=genotype),
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8),
             size=2, shape=21, stroke = 0.2, # stroke control border thickness
             color="black") +
  scale_fill_manual(values = rep("white", 2)) +
  ggbreak::scale_y_break(c(300, 2000)) +# check https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html#introduction
  annotate("text", x = c("Ebf1"), y = 2400, 
           label = paste(capture.output(as.data.frame(key_TF_genes_statistics_data_step7.2_3)), collapse = "\n"), # check https://stackoverflow.com/questions/14326573/print-a-data-frame-in-the-white-space-of-a-plot
           size=1.5, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('Transcription regulatory gene expression ( original TPM) \nin Trim66 KO by ggplot') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Transcription regulatory genes") + ylab("TPM (mean+-sem)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# save plot_step7.2_3 and plot_step7.2_4
pdf(paste(baseDir, "/results/key_TF_genes_TPM_Trim66_KO.pdf", sep = ""), paper = "a4")
plot_step7.2_3
plot_step7.2_4
dev.off()
# if you arrange plot_step7.2_3 and plot_step7.2_4 using ggarange, the y-break effect will be gone for unknown reason.
# You can also save plot_step7.2_3 and plot_step7.2_4 separately using ggsave, like this:
# ggsave(paste(baseDir, "/results/key_TF_genes_TPM_Trim66_KO_ggpubr::ggbarplot.pdf", sep=''), plot_step7.2_3, width = 25, height = 25, units ="cm")
# ggsave(paste(baseDir, "/results/key_TF_genes_TPM_Trim66_KO_ggplot.pdf", sep=''), plot_step7.2_4, width = 25, height = 25, units ="cm")

# 7.3 (1)TAAR genes using original TPM with ggpubr::ggbarplot
TAAR_genes_TPM_data_step7.3 <- TPM_final_for_save_step5[grepl("^Taar[2-9]", rownames(TPM_final_for_save_step5)), 1:10]
# statstical data of gene changes
TAAR_genes_statistics_data_step7.3 <- res_all[grepl("^Taar[2-9]", rownames(res_all)), c(1:2,5:6)]
TAAR_genes_statistics_data_step7.3$TPM_mean_het_5samples <- rowMeans(TAAR_genes_TPM_data_step7.3[,1:5])
TAAR_genes_statistics_data_step7.3$TPM_mean_homo_5samples <- rowMeans(TAAR_genes_TPM_data_step7.3[,6:10])
# arrange data for plotting
TAAR_genes_TPM_data_step7.3$name <- rownames(TAAR_genes_TPM_data_step7.3)
TAAR_genes_TPM_data_step7.3 <- melt(TAAR_genes_TPM_data_step7.3, id.vars = "name") #melt function from reshape2 package to convert dataframe to the format easier for ggplot
TAAR_genes_TPM_data_step7.3$genotype <- "het"
TAAR_genes_TPM_data_step7.3[grep("homo", TAAR_genes_TPM_data_step7.3$variable), ]$genotype <- "homo"
TAAR_genes_TPM_data_step7.3$genotype <- factor(TAAR_genes_TPM_data_step7.3$genotype, levels=c("het", "homo")) # force ggplot to plot het first
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
TAAR_genes_TPM_data_step7.3$name <- factor(TAAR_genes_TPM_data_step7.3$name, levels = rownames(TAAR_genes_statistics_data_step7.3))
#### use a different method to add text here, as using the above step 7.2 (1) method does not work for broken y axis with ggbreak::scale_y_break function
plot_step7.3_1 <- ggpubr::ggbarplot(TAAR_genes_TPM_data_step7.3, x = "name", y = "value", add = c("mean_se"),
                                    fill = "gray95", size = 0.25, # size to change the size of points and outlines. 
                                    add.params = list(size=.25), # Add.params to change the size for the argument "add", the error bar and point in this case
                                    title = "TAAR gene expression (original TPM) in Trim66 KO by ggpubr::ggbarplot", x.text.angle = 60, # rotate the x label
                                    ylab = "TPM (mean+-sem)", xlab = "Taar genes",
                                    color = "genotype", palette = c("springgreen3", "red1"),
                                    position = position_dodge(0.85)) + # to avoid add("mean_se", "jitter"), check https://github.com/kassambara/ggpubr/issues/130
  geom_point(aes(fill=genotype),
             size=2, shape=21, stroke = 0.2, # stroke control border thickness
             color="black",
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8)) +
  scale_fill_manual(values = rep("white", 2)) +
  annotate("text", x = c("Taar7d"), y = 100, 
           label = paste(capture.output(as.data.frame(TAAR_genes_statistics_data_step7.3)), collapse = "\n"), # check https://stackoverflow.com/questions/14326573/print-a-data-frame-in-the-white-space-of-a-plot
           size=2, colour = "black") # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console

#####################################################################################
#note: using ggplot2 is better
#####################################################################################
# 7.3 (2)TAAR genes using original TPM with ggplot
#### This is better. Next time use ggplot2 to make this kind of figures.
TAAR_TPM_means_se_step7.3_2 <- TAAR_genes_TPM_data_step7.3 %>%  group_by(name, genotype) %>% 
  summarize(mean_TPM=mean(value), 
            sd_TPM=sd(value), 
            Numbers=n(), 
            se=sd_TPM/sqrt(Numbers), 
            upper_limit=mean_TPM+se, 
            lower_limit=mean_TPM-se) 
TAAR_TPM_means_se_step7.3_2
# change name from character to factor according to the original gene order. Otherwise ggplot will change the order alphabetically.
TAAR_TPM_means_se_step7.3_2$name <- factor(TAAR_TPM_means_se_step7.3_2$name, levels = rownames(TAAR_genes_statistics_data_step7.3))
plot_step7.3_2 <- ggplot() + 
  geom_bar(TAAR_TPM_means_se_step7.3_2,
           mapping=aes(x=name, y=mean_TPM, color=genotype),
           stat="identity", 
           position = position_dodge(0.8),
           width = 0.7, #adjust the bar width
           size=0.25, # adjust the line thickness
           fill="gray95") + 
  geom_errorbar(TAAR_TPM_means_se_step7.3_2,
                mapping=aes(x=name, ymin=lower_limit, ymax=upper_limit, color=genotype),
                position = position_dodge(0.8), width = .2, size=0.25) +
  scale_color_manual(values = c("springgreen3", "red1")) +
  geom_point(TAAR_genes_TPM_data_step7.3,
             mapping=aes(x=name, y=value, group=genotype, fill=genotype),
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.8),
             size=2, shape=21, stroke = 0.2, # stroke control border thickness
             color="black") +
  scale_fill_manual(values = rep("white", 2)) +
  annotate("text", x = c("Taar7d"), y = 100, 
           label = paste(capture.output(as.data.frame(TAAR_genes_statistics_data_step7.3)), collapse = "\n"), # check https://stackoverflow.com/questions/14326573/print-a-data-frame-in-the-white-space-of-a-plot
           size=2, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('TAAR gene expression (original TPM) in Trim66 KO by ggplot') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("TAAR genes") + ylab("TPM (mean +- se)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

arrange_plots_step7.3 <- ggarrange(plot_step7.3_1, plot_step7.3_2, ncol = 1, nrow = 2)
ggsave(paste(baseDir, "/results/TAAR_genes_TPM_Trim66_KO.pdf", sep=''), arrange_plots_step7.3, width = 25, height = 25, units ="cm")

# plot percent TPM of all functional TAARs in het and homo
functional_TAAR_genes_TPM_data_step7.3_3 <- TAAR_genes_TPM_data_step7.3 %>% 
  filter(!grepl("Taar7c", name)) %>% 
  group_by(name, genotype) %>% 
  summarise(mean_TPM=mean(value)) %>% 
  tidyr::spread(key = genotype, value = mean_TPM) %>% #spread to convert long list to wide. check https://rpubs.com/mm-c/gather-and-spread
  rename(Trim66_KO_het_mean_TPM=het, Trim66_KO_homo_mean_TPM=homo) %>% 
  arrange(Trim66_KO_homo_mean_TPM) %>% 
  ungroup() %>% #ungroup the previous group_by name
  mutate(Trim66_KO_het_mean_TPM_percentage = Trim66_KO_het_mean_TPM/sum(Trim66_KO_het_mean_TPM), Trim66_KO_homo_mean_TPM_percentage = Trim66_KO_homo_mean_TPM/sum(Trim66_KO_homo_mean_TPM))
write.csv(functional_TAAR_genes_TPM_data_step7.3_3, file = paste(baseDir, "/results/Trim66_KO_het_vs_homo_RNA-seq_mean_TPM_functional_TAARs.csv", sep=''))
# what are the order of TAARs with descending TPM in homo
TAAR_order_TPM_desending_in_homo <- functional_TAAR_genes_TPM_data_step7.3_3 %>% 
  arrange(desc(Trim66_KO_homo_mean_TPM_percentage)) %>% 
  select(name) %>% 
  pull() %>% #pull to convert to vector
  paste(collapse = " ") # collapse to string
# what are the percentage of top 2 TAARs: TAAR2 and TAAR6
functional_TAAR_genes_TPM_data_step7.3_3 %>% 
  arrange(desc(Trim66_KO_homo_mean_TPM_percentage)) %>% 
  top_n(2) %>% 
  select(Trim66_KO_het_mean_TPM_percentage, Trim66_KO_homo_mean_TPM_percentage) %>% 
  colSums()
plot_step7.3_3 <-ggplot() + 
  geom_bar(functional_TAAR_genes_TPM_data_step7.3_3,
           mapping=aes(x="Trim66_KO_het", y=Trim66_KO_het_mean_TPM_percentage),
           stat="identity",
           color="black",
           fill = "white",
           width = 0.5, #adjust the bar width
           size=0.1) + # adjust the line thickness
  geom_bar(functional_TAAR_genes_TPM_data_step7.3_3,
           mapping=aes(x="Trim66_KO_homo", y=Trim66_KO_homo_mean_TPM_percentage),
           stat="identity",
           color="black",
           fill = "white",
           width = 0.5, #adjust the bar width
           size=0.1) + # adjust the line thickness
  theme(legend.position = "top",
        plot.title = element_text(size = rel(0.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1))) +
  ggtitle(paste('Percent TPM of all functional TAARs in het and homo \nsort by percent in homo, the order of TAARs is ', TAAR_order_TPM_desending_in_homo, '\ntop 2 TAARs (TAAR2 and TAAR6) occupy 83.50% in homo and 24.69% in het', sep = "")) +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM percentage") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/functional_TAAR_genes_TPM_percentage_Trim66_KO.pdf", sep=''), plot_step7.3_3, width = 25, height = 25, units ="cm")

# 7.4 (1)functional OR genes (all together) using original TPM
functional_OR_genes_TPM_data_step7.4_1 <- TPM_final_for_save_step5[grepl("^Olfr", rownames(TPM_final_for_save_step5)) & grepl("protein_coding", res_all$gene_type), 1:10]
functional_OR_genes_TPM_data_step7.4_1$name <- rownames(functional_OR_genes_TPM_data_step7.4_1)
functional_OR_genes_TPM_data_step7.4_1 <- melt(functional_OR_genes_TPM_data_step7.4_1, id.vars = "name") #melt function from reshape2 package to convert dataframe to the format easier for ggplot
functional_OR_genes_TPM_data_step7.4_1$genotype <- "het"
functional_OR_genes_TPM_data_step7.4_1[grep("homo", functional_OR_genes_TPM_data_step7.4_1$variable), ]$genotype <- "homo"
functional_OR_genes_TPM_data_step7.4_1$genotype <- factor(functional_OR_genes_TPM_data_step7.4_1$genotype, levels=c("het", "homo")) # force ggplot to plot het first
# unpaired t.test
t.test_unpaired_functional_OR_genes_step7.4_1 <- t.test(functional_OR_genes_TPM_data_step7.4_1[functional_OR_genes_TPM_data_step7.4_1$genotype =="homo", "value"], 
                                                        functional_OR_genes_TPM_data_step7.4_1[functional_OR_genes_TPM_data_step7.4_1$genotype =="het", "value"])
# This plot with mean +- se is not ideal, as the mean and se are very low
functional_OR_genes_TPM_means_se_step7.4_1 <- functional_OR_genes_TPM_data_step7.4_1 %>%  group_by(genotype) %>% 
  summarize(mean_TPM=mean(value), 
            sd_TPM=sd(value), 
            Numbers=n(), 
            se=sd_TPM/sqrt(Numbers), 
            upper_limit=mean_TPM+se, 
            lower_limit=mean_TPM-se,
            sum_TPM=sum(value)) 
functional_OR_genes_TPM_means_se_step7.4_1
plot_step7.4_1_1 <-ggplot() + 
  geom_point(functional_OR_genes_TPM_data_step7.4_1,
             mapping=aes(x=genotype, y=value),
             position=position_jitter(width = 0.2),
             size=1, shape=21, alpha=0.4, 
             stroke = 0.1, # stroke control border thickness
             color="gray40", fill=NA) +
  geom_bar(functional_OR_genes_TPM_means_se_step7.4_1,
           mapping=aes(x=genotype, y=mean_TPM, color=genotype),
           stat="identity", 
           width = 0.7, #adjust the bar width
           size=0.1, # adjust the line thickness
           fill=NA) + 
  geom_errorbar(functional_OR_genes_TPM_means_se_step7.4_1,
                mapping=aes(x=genotype, ymin=lower_limit, ymax=upper_limit, color=genotype),
                width = .5, size=0.1) +
  scale_color_manual(values = c("springgreen3", "red1")) +
  annotate("text", x = c("het"), y = 400, 
           label = paste("unpaired Student's t-test, p value: ", 
                         t.test_unpaired_functional_OR_genes_step7.4_1$p.value, sep = ""),
           size=2, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('Functional OR gene expression (original TPM, mean+-se) in Trim66 KO') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM (mean +- se)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/functional_OR_genes_TPM_Trim66_KO.pdf", sep=''), plot_step7.4_1_1, width = 25, height = 25, units ="cm")

# use ggplot to plot the median and the 25th and 75th percentiles and the upper and lower whiskers
plot_step7.4_1_2 <- ggplot(functional_OR_genes_TPM_data_step7.4_1, aes(x = genotype, y = value)) + 
  geom_boxplot(aes(fill = genotype), outlier.shape=NA) + # check https://stackoverflow.com/questions/34602872/how-to-remove-dots-and-extend-boxplots-in-ggplot2
  scale_fill_manual(values = c("springgreen3", "lightcoral")) +
  ggbreak::scale_y_break(c(15, 450), scales = 0.01) +#use ggbreak::scale_y_break instead of set ylim, as set ylim will remove several points.  check https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html#introduction
  annotate("text", x = c("homo"), y = 10, 
           label = paste("unpaired Student's t-test, p value: ", t.test_unpaired_functional_OR_genes_step7.4_1$p.value, sep = ""),
           size=10, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional OR gene expression (original TPM, median+-quantile) in Trim66 KO by ggplot with y break') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM (median+-the first and third quartiles (the 25th and 75th percentiles) \nThe upper whisker extends from the hinge to the largest value no further than 1.5 * IQR from the hinge \n(where IQR is the inter-quartile range, or distance between the first and third quartiles). \nThe lower whisker extends from the hinge to the smallest value at most 1.5 * IQR of the hinge.)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

ggsave(paste(baseDir, "/results/functional_OR_genes_TPM_Trim66_KO_ggpplot_with_y_break.pdf", sep=''), plot_step7.4_1_2, width = 100, height = 100, units ="cm")

# sum of TPM from all functional ORs
plot_step7.4_1_3 <-ggplot() + 
  geom_bar(functional_OR_genes_TPM_means_se_step7.4_1,
           mapping=aes(x=genotype, y=sum_TPM, color=genotype, fill=genotype),
           stat="identity",
           fill = "gray95",
           width = 0.7, #adjust the bar width
           size=0.25) + # adjust the line thickness
  scale_color_manual(values = c("springgreen3", "red1")) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ylim(0,40000) +
  ggtitle('Functional OR gene expression (original TPM, sum) in Trim66 KO') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM (all functional ORs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/functional_OR_genes_TPM_sum_Trim66_KO.pdf", sep=''), plot_step7.4_1_3, width = 25, height = 25, units ="cm")

# plot percent TPM of all functional ORs in het and homo
functional_OR_genes_TPM_data_step7.4_1_4 <- functional_OR_genes_TPM_data_step7.4_1 %>% 
  group_by(name, genotype) %>% 
  summarise(mean_TPM=mean(value)) %>% 
  tidyr::spread(key = genotype, value = mean_TPM) %>% #spread to convert long list to wide. check https://rpubs.com/mm-c/gather-and-spread
  rename(Trim66_KO_het_mean_TPM=het, Trim66_KO_homo_mean_TPM=homo) %>% 
  arrange(Trim66_KO_homo_mean_TPM) %>% 
  ungroup() %>% #ungroup the previous group_by name
  mutate(Trim66_KO_het_mean_TPM_percentage = Trim66_KO_het_mean_TPM/sum(Trim66_KO_het_mean_TPM), Trim66_KO_homo_mean_TPM_percentage = Trim66_KO_homo_mean_TPM/sum(Trim66_KO_homo_mean_TPM))
write.csv(functional_OR_genes_TPM_data_step7.4_1_4, file = paste(baseDir, "/results/Trim66_KO_het_vs_homo_RNA-seq_mean_TPM_functional_ORs.csv", sep=''))
# what are the top 10 ORs with most TPM in homo
top_10_ORs_TPM_in_homo <- functional_OR_genes_TPM_data_step7.4_1_4 %>% 
  arrange(desc(Trim66_KO_homo_mean_TPM_percentage)) %>% 
  top_n(10) %>% 
  select(name) %>% 
  pull() %>% #pull to convert to vector
  paste(collapse = " ") # collapse to string
# what are the percentage of top 10 ORs
functional_OR_genes_TPM_data_step7.4_1_4 %>% 
  arrange(desc(Trim66_KO_homo_mean_TPM_percentage)) %>% 
  top_n(10) %>% 
  select(Trim66_KO_het_mean_TPM_percentage, Trim66_KO_homo_mean_TPM_percentage) %>% 
  colSums()
plot_step7.4_1_4 <-ggplot() + 
  geom_bar(functional_OR_genes_TPM_data_step7.4_1_4,
           mapping=aes(x="Trim66_KO_het", y=Trim66_KO_het_mean_TPM_percentage),
           stat="identity",
           color="black",
           fill = "white",
           width = 0.5, #adjust the bar width
           size=0.1) + # adjust the line thickness
  geom_bar(functional_OR_genes_TPM_data_step7.4_1_4,
           mapping=aes(x="Trim66_KO_homo", y=Trim66_KO_homo_mean_TPM_percentage),
           stat="identity",
           color="black",
           fill = "white",
           width = 0.5, #adjust the bar width
           size=0.1) + # adjust the line thickness
  theme(legend.position = "top",
        plot.title = element_text(size = rel(0.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1))) +
  ggtitle(paste('Percent TPM of all functional ORs in het and homo \nsort by percent in homo, top 10 ORs are ', top_10_ORs_TPM_in_homo, '\ntop 10 ORs occupy 47.97% in homo and 3.13% in het', sep = "")) +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM percentage") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/functional_OR_genes_TPM_percentage_Trim66_KO.pdf", sep=''), plot_step7.4_1_4, width = 25, height = 25, units ="cm")

# 7.4 (2)pseudogene OR genes (all together) using original TPM
pseudogene_OR_genes_TPM_data_step7.4_2 <- TPM_final_for_save_step5[grepl("^Olfr", rownames(TPM_final_for_save_step5)) & grepl("pseudogene", res_all$gene_type), 1:10]
pseudogene_OR_genes_TPM_data_step7.4_2$name <- rownames(pseudogene_OR_genes_TPM_data_step7.4_2)
pseudogene_OR_genes_TPM_data_step7.4_2 <- melt(pseudogene_OR_genes_TPM_data_step7.4_2, id.vars = "name") #melt function from reshape2 package to convert dataframe to the format easier for ggplot
pseudogene_OR_genes_TPM_data_step7.4_2$genotype <- "het"
pseudogene_OR_genes_TPM_data_step7.4_2[grep("homo", pseudogene_OR_genes_TPM_data_step7.4_2$variable), ]$genotype <- "homo"
pseudogene_OR_genes_TPM_data_step7.4_2$genotype <- factor(pseudogene_OR_genes_TPM_data_step7.4_2$genotype, levels=c("het", "homo")) # force ggplot to plot het first
# unpaired t.test
t.test_unpaired_pseudogene_OR_genes_step7.4_2 <- t.test(pseudogene_OR_genes_TPM_data_step7.4_2[pseudogene_OR_genes_TPM_data_step7.4_2$genotype =="homo", "value"], 
                                                        pseudogene_OR_genes_TPM_data_step7.4_2[pseudogene_OR_genes_TPM_data_step7.4_2$genotype =="het", "value"])
# This plot with mean +- se is not ideal, as the mean and se are very low
pseudogene_OR_genes_TPM_means_se_step7.4_2 <- pseudogene_OR_genes_TPM_data_step7.4_2 %>%  group_by(genotype) %>% 
  summarize(mean_TPM=mean(value), 
            sd_TPM=sd(value), 
            Numbers=n(), 
            se=sd_TPM/sqrt(Numbers), 
            upper_limit=mean_TPM+se, 
            lower_limit=mean_TPM-se,
            sum_TPM=sum(value)) 
pseudogene_OR_genes_TPM_means_se_step7.4_2
plot_step7.4_2_1 <-ggplot() + 
  geom_point(pseudogene_OR_genes_TPM_data_step7.4_2,
             mapping=aes(x=genotype, y=value),
             position=position_jitter(width = 0.2),
             size=1, shape=21, alpha=0.4, 
             stroke = 0.1, # stroke control border thickness
             color="gray40", fill=NA) +
  geom_bar(pseudogene_OR_genes_TPM_means_se_step7.4_2,
           mapping=aes(x=genotype, y=mean_TPM, color=genotype),
           stat="identity", 
           width = 0.7, #adjust the bar width
           size=0.1, # adjust the line thickness
           fill=NA) + 
  geom_errorbar(pseudogene_OR_genes_TPM_means_se_step7.4_2,
                mapping=aes(x=genotype, ymin=lower_limit, ymax=upper_limit, color=genotype),
                width = .5, size=0.1) +
  scale_color_manual(values = c("springgreen3", "red1")) +
  annotate("text", x = c("het"), y = 180, 
           label = paste("unpaired Student's t-test, p value: ", 
                         round(t.test_unpaired_pseudogene_OR_genes_step7.4_2$p.value, 4), sep = ""),
           size=2, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('pseudogene OR gene expression (original TPM, mean+-se) in Trim66 KO') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM (mean +- se)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/pseudogene_OR_genes_TPM_Trim66_KO.pdf", sep=''), plot_step7.4_2_1, width = 25, height = 25, units ="cm")

# use ggplot to plot the median and the 25th and 75th percentiles and the upper and lower whiskers
plot_step7.4_2_2 <- ggplot(pseudogene_OR_genes_TPM_data_step7.4_2, aes(x = genotype, y = value)) + 
  geom_boxplot(aes(fill = genotype), outlier.shape=NA) + # check https://stackoverflow.com/questions/34602872/how-to-remove-dots-and-extend-boxplots-in-ggplot2
  scale_fill_manual(values = c("springgreen3", "lightcoral")) +
  ggbreak::scale_y_break(c(1, 200), scales = 0.01) +#use ggbreak::scale_y_break instead of set ylim, as set ylim will remove several points.  check https://cran.r-project.org/web/packages/ggbreak/vignettes/ggbreak.html#introduction
  annotate("text", x = c("homo"), y = 0.75, 
           label = paste("unpaired Student's t-test, p value: ", round(t.test_unpaired_pseudogene_OR_genes_step7.4_2$p.value,4), sep = ""),
           size=5, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('pseudogene OR gene expression (original TPM, median+-quantile) in Trim66 KO by ggplot with y break') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM (median+-the first and third quartiles (the 25th and 75th percentiles) \nThe upper whisker extends from the hinge to the largest value no further than 1.5 * IQR from the hinge \n(where IQR is the inter-quartile range, or distance between the first and third quartiles). \nThe lower whisker extends from the hinge to the smallest value at most 1.5 * IQR of the hinge.)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/pseudogene_OR_genes_TPM_Trim66_KO_ggpplot_with_y_break.pdf", sep=''), plot_step7.4_2_2, width = 25, height = 25, units ="cm")

# sum of TPM from all pseudogene OR genes
plot_step7.4_2_3 <-ggplot() + 
  geom_bar(pseudogene_OR_genes_TPM_means_se_step7.4_2,
           mapping=aes(x=genotype, y=sum_TPM, color=genotype, fill=genotype),
           stat="identity", 
           width = 0.7, #adjust the bar width
           size=0.1) + # adjust the line thickness
  scale_color_manual(values = c("springgreen3", "red1")) +
  scale_fill_manual(values = c("springgreen3", "red1")) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(1.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1)),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  ggtitle('pseudogene OR gene expression (original TPM, sum) in Trim66 KO') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM (all pseudogene ORs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/pseudogene_OR_genes_TPM_sum_Trim66_KO.pdf", sep=''), plot_step7.4_2_3, width = 25, height = 25, units ="cm")

# plot percent TPM of all pseudogene ORs in het and homo
pseudogene_OR_genes_TPM_data_step7.4_2_4 <- pseudogene_OR_genes_TPM_data_step7.4_2 %>% 
  group_by(name, genotype) %>% 
  summarise(mean_TPM=mean(value)) %>% 
  tidyr::spread(key = genotype, value = mean_TPM) %>% #spread to convert long list to wide. check https://rpubs.com/mm-c/gather-and-spread
  rename(Trim66_KO_het_mean_TPM=het, Trim66_KO_homo_mean_TPM=homo) %>% 
  arrange(Trim66_KO_homo_mean_TPM) %>% 
  ungroup() %>% #ungroup the previous group_by name
  mutate(Trim66_KO_het_mean_TPM_percentage = Trim66_KO_het_mean_TPM/sum(Trim66_KO_het_mean_TPM), Trim66_KO_homo_mean_TPM_percentage = Trim66_KO_homo_mean_TPM/sum(Trim66_KO_homo_mean_TPM))
write.csv(pseudogene_OR_genes_TPM_data_step7.4_2_4, file = paste(baseDir, "/results/Trim66_KO_het_vs_homo_RNA-seq_mean_TPM_pesudogene_ORs.csv", sep=''))
# what are the top 10 ORs with most TPM in homo
top_10_ORs_TPM_in_homo <- pseudogene_OR_genes_TPM_data_step7.4_2_4 %>% 
  arrange(desc(Trim66_KO_homo_mean_TPM_percentage)) %>% 
  top_n(10) %>% 
  select(name) %>% 
  pull() %>% #pull to convert to vector
  paste(collapse = " ") # collapse to string
# what are the percentage of top 10 ORs
pseudogene_OR_genes_TPM_data_step7.4_2_4 %>% 
  arrange(desc(Trim66_KO_homo_mean_TPM_percentage)) %>% 
  top_n(10) %>% 
  select(Trim66_KO_het_mean_TPM_percentage, Trim66_KO_homo_mean_TPM_percentage) %>% 
  colSums()
plot_step7.4_1_4 <-ggplot() + 
  geom_bar(pseudogene_OR_genes_TPM_data_step7.4_2_4,
           mapping=aes(x="Trim66_KO_het", y=Trim66_KO_het_mean_TPM_percentage),
           stat="identity",
           color="black",
           fill = "white",
           width = 0.5, #adjust the bar width
           size=0.1) + # adjust the line thickness
  geom_bar(pseudogene_OR_genes_TPM_data_step7.4_2_4,
           mapping=aes(x="Trim66_KO_homo", y=Trim66_KO_homo_mean_TPM_percentage),
           stat="identity",
           color="black",
           fill = "white",
           width = 0.5, #adjust the bar width
           size=0.1) + # adjust the line thickness
  theme(legend.position = "top",
        plot.title = element_text(size = rel(0.5)), # adjust the title size
        axis.title = element_text(size = rel(1)),
        axis.text = element_text(size = rel(1))) +
  ggtitle(paste('Percent TPM of all pseudogene ORs in het and homo \nsort by percent in homo, top 10 ORs are ', top_10_ORs_TPM_in_homo, '\ntop 10 ORs occupy 82.17% in homo and 34.06% in het', sep = "")) +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Genotype") + ylab("TPM percentage") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/
ggsave(paste(baseDir, "/results/pseudogene_OR_genes_TPM_percentage_Trim66_KO.pdf", sep=''), plot_step7.4_1_4, width = 25, height = 25, units ="cm")

#####################################################################################################################
       # Summary of best way for plotting points and mean+-se or boxplot or violin
# 1. to plot points and mean+-se together, firstly generate a new data frame with mean and se information using dplyr::summarize
      # then use geom_point to plot points, geom_bar to plot mean, geom_errorbar to plot errorbar (see above code).
# 2. to plot points and boxplot (meadian+-quantile) or violin, do not need to generate new dataframe with mean and se.
      # instead, use geom_boxplot or geom_violin to plot boxplot or violin, then use geom_point to plot points on top.
# 3. note that put geom_point as last if you want points on top layer. On the contrary, put geom_point first if you want points at bottom layer.
#####################################################################################################################



# Step 8: Circos plot for OR genes changes --------------------------------


#####################################################################################################################
      # Plot Step 8: plot Circos plot for OR gene changes
      # https://jokergoo.github.io/circlize_book/book/introduction.html
      # https://cran.r-project.org/web/packages/circlize/circlize.pdf
#####################################################################################################################
# install.packages("circlize")
library(circlize)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap) # to add legend in Circos plot

# step 8.1: Circos plot of functional OR genes --------------------------------------

#8.1 plot 1: het (green bar) and homo (red bar) use a single file with combined het and homo mean TPM
# here use het as green and homo as red to keep consistent with the above rule: downregulation as green and upregulation as red
# extract the functional OR data
TPM_final_for_save_step8.1 <- TPM_final_for_save_step5
functional_OR_genes_TPM_data_for_Circos_step8.1 <- TPM_final_for_save_step8.1[grepl("^Olfr", rownames(TPM_final_for_save_step8.1)) & grepl("protein_coding", res_all$gene_type), 1:10]
functional_OR_genes_TPM_data_for_Circos_step8.1$het_TPM_5samples <- rowMeans(functional_OR_genes_TPM_data_for_Circos_step8.1[,1:5])
functional_OR_genes_TPM_data_for_Circos_step8.1$homo_TPM_5samples <- rowMeans(functional_OR_genes_TPM_data_for_Circos_step8.1[,6:10])
# obtain het data and make relative chr coordinates (start from 1 to sum, end = start+0.8) for plotting with gap 0.2
# the order of column should be similar to bed file, i.e., chr name, chr_start, chr_end, value1, value2....
functional_OR_genes_TPM_data_for_Circos_step8.1_1 <- res_functional_OR["chr_name"]
functional_OR_genes_TPM_data_for_Circos_step8.1_1$genomic_start_relative <- NA
for (i in unique(functional_OR_genes_TPM_data_for_Circos_step8.1_1$chr_name)) {
  functional_OR_genes_TPM_data_for_Circos_step8.1_1[functional_OR_genes_TPM_data_for_Circos_step8.1_1$chr_name == i, ]$genomic_start_relative <- 1:sum(functional_OR_genes_TPM_data_for_Circos_step8.1_1$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}
functional_OR_genes_TPM_data_for_Circos_step8.1_1$genomic_end_relative <- functional_OR_genes_TPM_data_for_Circos_step8.1_1$genomic_start_relative + 0.8
functional_OR_genes_TPM_data_for_Circos_step8.1_1$TPM_5samples <- functional_OR_genes_TPM_data_for_Circos_step8.1$het_TPM_5samples
functional_OR_genes_TPM_data_for_Circos_step8.1_1$genotype <- factor("het")
# same to obtain homo data
functional_OR_genes_TPM_data_for_Circos_step8.1_2 <- res_functional_OR["chr_name"]
functional_OR_genes_TPM_data_for_Circos_step8.1_2$genomic_start_relative <- NA
for (i in unique(functional_OR_genes_TPM_data_for_Circos_step8.1_2$chr_name)) {
  functional_OR_genes_TPM_data_for_Circos_step8.1_2[functional_OR_genes_TPM_data_for_Circos_step8.1_2$chr_name == i, ]$genomic_start_relative <- 1:sum(functional_OR_genes_TPM_data_for_Circos_step8.1_2$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}
functional_OR_genes_TPM_data_for_Circos_step8.1_2$genomic_end_relative <- functional_OR_genes_TPM_data_for_Circos_step8.1_2$genomic_start_relative + 0.8
functional_OR_genes_TPM_data_for_Circos_step8.1_2$TPM_5samples <- functional_OR_genes_TPM_data_for_Circos_step8.1$homo_TPM_5samples
functional_OR_genes_TPM_data_for_Circos_step8.1_2$genotype <- factor("homo")
# rbind the het and homo data
functional_OR_genes_TPM_data_for_Circos_final_step8.1 <- rbind(functional_OR_genes_TPM_data_for_Circos_step8.1_1, functional_OR_genes_TPM_data_for_Circos_step8.1_2)
functional_OR_genes_TPM_data_for_Circos_final_step8.1 <- as.data.frame(functional_OR_genes_TPM_data_for_Circos_final_step8.1)

# not using scale_color_brewer as there are 21 different chr, but most color palettes have 8-12 colors. check http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html  and https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_chr_step8.1 <- col_for_chr_universal
# obtain the colors for chr used in step8.1
col_for_chr_step8.1 <- col_for_chr_step8.1[unique(functional_OR_genes_TPM_data_for_Circos_final_step8.1$chr_name)]
# start plotting plot 1
pdf(paste(baseDir, "/results/functional_OR_genes_TPM_Trim66_KO_Circos_plot.pdf", sep=''), paper='a4')
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1, plotType = "labels")
circos.initializeWithIdeogram(functional_OR_genes_TPM_data_for_Circos_final_step8.1,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.1,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(functional_OR_genes_TPM_data_for_Circos_final_step8.1,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border=NA, #otherwise there are black border
                                         col = ifelse(value[[2]] =="het", "springgreen3", "lightcoral")) #use [[]] instead of [] to obtain vector other than dataframe.
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) \nof functional OR genes in Trim66 KO",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "original TPM value", cex = 0.75)
# add legend
# check Circos mannual (https://jokergoo.github.io/circlize_book/book/index.html) and ComplexHeatmap mannual (https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html)
legend_boxs_step8.1 <- ComplexHeatmap::Legend(labels = c("het", "homo"), 
                                             type = "boxplot", 
                                             legend_gp = gpar(fill=c("springgreen3", "lightcoral")), 
                                             title_position = "topleft", 
                                             title = "Genotype")
# "npc" means Normalised Parent Coordinates (the default for unit). The origin of the viewport is (0, 0) and the viewport has a width and height of 1 unit. For example, (0.5, 0.5) is the centre of the viewport.
# check https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/unit
ComplexHeatmap::draw(legend_boxs_step8.1, x = unit(0.1, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
circos.clear()

#8.1 plot 2: het (green bar) and homo (red bar) with TPM above 50 changed to 50 for better visualization
quantile(functional_OR_genes_TPM_data_for_Circos_final_step8.1$TPM_5samples, c(0.5, 0.75, 0.9,0.95, 0.98,0.99,1))
# since 98%, and 99% percentil of the TPM values is 27.99 and 42.15, I will manipulate the TPM value above 50 to be 50 for better visualizaiton.
sum(functional_OR_genes_TPM_data_for_Circos_final_step8.1$TPM_5samples>50)  #total 18 data points with TPM above 50
functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1 <- functional_OR_genes_TPM_data_for_Circos_final_step8.1
functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1[functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1$TPM_5samples>50, ]$TPM_5samples <- 50
# start plotting plot 2
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1, plotType = "labels")
circos.initializeWithIdeogram(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.1,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border=NA, #otherwise there are black border
                                         col = ifelse(value[[2]] =="het", "springgreen3", "lightcoral")) #use [[]] instead of [] to obtain vector other than dataframe.
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) \nof functional OR genes in Trim66 KO",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "TPM larger than 50 are \nchanged to 50 for better visualization\n(18 data points out of 1137*2 are changed)", cex = 0.75)
# add legend
ComplexHeatmap::draw(legend_boxs_step8.1, x = unit(0.1, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
circos.clear()
dev.off()


# step 8.2: Circos plot of functional OR genes with significant changes marked --------

#8.2: plot 1- het_significant (green bar), homo_significant (red bar), het_nonsig (gray 80 bar) and homo_nonsig (gray10 bar) use a single file with combined het and homo mean TPM
# extract the functional OR data
TPM_final_for_save_step8.2 <- TPM_final_for_save_step8.1
functional_OR_genes_TPM_data_for_Circos_step8.2 <- TPM_final_for_save_step8.2[grepl("^Olfr", rownames(TPM_final_for_save_step8.2)) & grepl("protein_coding", res_all$gene_type), 1:10]
functional_OR_genes_TPM_data_for_Circos_step8.2$het_TPM_5samples <- rowMeans(functional_OR_genes_TPM_data_for_Circos_step8.2[,1:5])
functional_OR_genes_TPM_data_for_Circos_step8.2$homo_TPM_5samples <- rowMeans(functional_OR_genes_TPM_data_for_Circos_step8.2[,6:10])
# obtain het data and make relative chr coordinates (start from 1 to sum, end = start+0.8)
functional_OR_genes_TPM_data_for_Circos_step8.2_1 <- res_functional_OR["chr_name"]
functional_OR_genes_TPM_data_for_Circos_step8.2_1$genomic_start_relative <- NA
for (i in unique(functional_OR_genes_TPM_data_for_Circos_step8.2_1$chr_name)) {
  functional_OR_genes_TPM_data_for_Circos_step8.2_1[functional_OR_genes_TPM_data_for_Circos_step8.2_1$chr_name == i, ]$genomic_start_relative <- 1:sum(functional_OR_genes_TPM_data_for_Circos_step8.2_1$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}
functional_OR_genes_TPM_data_for_Circos_step8.2_1$genomic_end_relative <- functional_OR_genes_TPM_data_for_Circos_step8.2_1$genomic_start_relative + 0.8
functional_OR_genes_TPM_data_for_Circos_step8.2_1$TPM_5samples <- functional_OR_genes_TPM_data_for_Circos_step8.2$het_TPM_5samples
functional_OR_genes_TPM_data_for_Circos_step8.2_1$genotype <- factor("het")
functional_OR_genes_TPM_data_for_Circos_step8.2_1$is_significant_padj_0.05_log2FC_0.585 <- res_functional_OR$padj < 0.05 & abs(res_functional_OR$log2FoldChange)>=0.585
functional_OR_genes_TPM_data_for_Circos_step8.2_1$genotype_with_is_significant <- paste(functional_OR_genes_TPM_data_for_Circos_step8.2_1$genotype, 
                                                                                        functional_OR_genes_TPM_data_for_Circos_step8.2_1$is_significant_padj_0.05_log2FC_0.585, 
                                                                                        sep = "_")
# same to obtain homo data
functional_OR_genes_TPM_data_for_Circos_step8.2_2 <- res_functional_OR["chr_name"]
functional_OR_genes_TPM_data_for_Circos_step8.2_2$genomic_start_relative <- NA
for (i in unique(functional_OR_genes_TPM_data_for_Circos_step8.2_2$chr_name)) {
  functional_OR_genes_TPM_data_for_Circos_step8.2_2[functional_OR_genes_TPM_data_for_Circos_step8.2_2$chr_name == i, ]$genomic_start_relative <- 1:sum(functional_OR_genes_TPM_data_for_Circos_step8.2_2$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}
functional_OR_genes_TPM_data_for_Circos_step8.2_2$genomic_end_relative <- functional_OR_genes_TPM_data_for_Circos_step8.2_2$genomic_start_relative + 0.8
functional_OR_genes_TPM_data_for_Circos_step8.2_2$TPM_5samples <- functional_OR_genes_TPM_data_for_Circos_step8.2$homo_TPM_5samples
functional_OR_genes_TPM_data_for_Circos_step8.2_2$genotype <- factor("homo")
functional_OR_genes_TPM_data_for_Circos_step8.2_2$is_significant_padj_0.05_log2FC_0.585 <- res_functional_OR$padj < 0.05 & abs(res_functional_OR$log2FoldChange)>=0.585
functional_OR_genes_TPM_data_for_Circos_step8.2_2$genotype_with_is_significant <- paste(functional_OR_genes_TPM_data_for_Circos_step8.2_2$genotype, 
                                                                                        functional_OR_genes_TPM_data_for_Circos_step8.2_2$is_significant_padj_0.05_log2FC_0.585, 
                                                                                        sep = "_")
# rbind the het and homo data
functional_OR_genes_TPM_data_for_Circos_final_step8.2 <- rbind(functional_OR_genes_TPM_data_for_Circos_step8.2_1, functional_OR_genes_TPM_data_for_Circos_step8.2_2)
functional_OR_genes_TPM_data_for_Circos_final_step8.2 <- as.data.frame(functional_OR_genes_TPM_data_for_Circos_final_step8.2)
# start plotting plot 1
pdf(paste(baseDir, "/results/functional_OR_genes_TPM_Trim66_KO_Circos_plot_sig_nonsig_labeled.pdf", sep=''), paper = 'a4')
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_step8.2, plotType = "labels")
circos.initializeWithIdeogram(functional_OR_genes_TPM_data_for_Circos_final_step8.2,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
col_for_chr_step8.2 <- col_for_chr_step8.1
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.2,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
library(hutils) # to use Switch, which is the vectorized form of base::switch. like ifelse is vevtorized form of if else. Check: https://www.delftstack.com/howto/r/use-a-vectorized-if-function-with-multiple-conditions-in-r/
circos.genomicTrack(functional_OR_genes_TPM_data_for_Circos_final_step8.2,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data.
                                         border=NA, #otherwise there are black border
                                         col = hutils::Switch(value[[4]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                              het_TRUE="springgreen3", 
                                                              homo_TRUE="lightcoral",
                                                              het_FALSE="gray80",
                                                              homo_FALSE="gray10",
                                                              DEFAULT = "white"))
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) of functional OR genes in Trim66 KO \nwith information of signicant changes padj<0.05 and abs(log2FC)>=0.585",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "original TPM value", cex = 0.75)
# add legend
# check Circos mannual (https://jokergoo.github.io/circlize_book/book/index.html) and ComplexHeatmap mannual (https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html)
legend_boxs_step8.2 <- ComplexHeatmap::Legend(labels = c("het sig", "homo sig", "het non_sig", "homo_non_sig"), 
                                             type = "boxplot", 
                                             legend_gp = gpar(fill=c("springgreen3", "lightcoral", "gray80", "gray10")), 
                                             title_position = "topleft", 
                                             title = "Genotype")
# "npc" means Normalised Parent Coordinates (the default for unit). The origin of the viewport is (0, 0) and the viewport has a width and height of 1 unit. For example, (0.5, 0.5) is the centre of the viewport.
# check https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/unit
ComplexHeatmap::draw(legend_boxs_step8.2, 
                     x = unit(0.05, "npc"), 
                     y = unit(0.92, "npc"), 
                     just = c("left", "top"))
circos.clear()

# 8.2 plot 2: since 98%, and 99% percentil of the TPM values is 27.99 and 42.14, I will manipulate the TPM value above 50 to be 50 for better visualizaiton.
sum(functional_OR_genes_TPM_data_for_Circos_final_step8.2$TPM_5samples>50)  #total 18 data points with TPM above 50
functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.2 <- functional_OR_genes_TPM_data_for_Circos_final_step8.2
functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.2[functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.2$TPM_5samples>50, ]$TPM_5samples <- 50
# start plotting plot 2
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.2, plotType = "labels")
circos.initializeWithIdeogram(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.2,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.2,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.2,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data.
                                         border=NA, #otherwise there are black border
                                         col = hutils::Switch(value[[4]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                              het_TRUE="springgreen3", 
                                                              homo_TRUE="lightcoral",
                                                              het_FALSE="gray80",
                                                              homo_FALSE="gray10",
                                                              DEFAULT = "white"))
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) of functional OR genes in Trim66 KO \nwith information of signicant changes padj<0.05 and abs(log2FC)>=0.585",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "TPM larger than 50 are \nchanged to 50 for better visualization\n(18 data points out of 1137*2 are changed)", cex = 0.75)
# add legend
ComplexHeatmap::draw(legend_boxs_step8.2, 
                     x = unit(0.05, "npc"), 
                     y = unit(0.92, "npc"), 
                     just = c("left", "top"))
circos.clear()
dev.off()



# !! Summary and tips for Circos plot (do not run) !!! --------------------
   Note: There are some better ways to plot in this summary!!!
#####################################################################################################
     # Summary and tips for Circos plot from example of functional OR gene dataset
     #  (read the online manual carefully)
     #    https://jokergoo.github.io/circlize_book/book/index.html
###################################################################################################### Aim: set the color for different chromosomes, which is consistent with plot step 4.1 color setup.
#####################################################################################################
# step 8.1 tips 1: to set up a set of colors and show the colors
   To show a set of customized color in the figure.
###################################################################################################### Aim: set the color for different chromosomes, which is consistent with plot step 4.1 color setup.
# not using scale_color_brewer as there are 21 different chr, but most color palettes have 8-12 colors. check http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html  and https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_chr_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21) ## not using scale_color_brewer as there are 21 different chr, but most color palettes have 8-12 colors. check http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html  and https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
names(col_for_chr_universal) <- unique(res_all$chr_name)
# to show the different colors with points
plot(NULL, xlim = c(1, length(col_for_chr_universal)), ylim = c(0, 2), axes = FALSE, ann = FALSE)
points(1:length(names(col_for_chr_universal)), rep(1, length(names(col_for_chr_universal))), 
       pch = 16, cex = 3, 
       col = col_for_chr_universal)
text(1:length(names(col_for_chr_universal)), rep(1.2, length(names(col_for_chr_universal))),
     names(col_for_chr_universal),
     cex=0.6)
#####################################################################################################

#####################################################################################################
# step 8.1 tips 2: to plot chromosomes with different colors
     The above code is better. This is just for record!
# For the record of different ways to plot chromosomes with different colors similar to above circos.track()
###################################################################################################### another way 1:
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.numeric.index # check CELL_META information. https://jokergoo.github.io/circlize_book/book/circular-layout.html
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  l = CELL_META$sector.index == names(col_for_chr_step8.1)
  circos.rect(xlim[1], 0, xlim[2], 1, col = col_for_chr_step8.1[l])
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.05, bg.border = NA)
# another way 2: directly obtain data using index of chr variable
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = col_for_chr_step8.1[chr])
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.05, bg.border = NA)
#####################################################################################################

####################################################################################
   # Step 8.1 tip3: Different way 1 to plot (functional OR as example): Use list file
           The above code is better. This is just for record!
####################################################################################
functional_OR_genes_TPM_data_for_Circos_step8.1_1 <- as.data.frame(functional_OR_genes_TPM_data_for_Circos_step8.1_1)
functional_OR_genes_TPM_data_for_Circos_step8.1_2 <- as.data.frame(functional_OR_genes_TPM_data_for_Circos_step8.1_2)
functional_OR_genes_TPM_data_for_Circos_step8.1_list = list(functional_OR_genes_TPM_data_for_Circos_step8.1_1,
                                                            functional_OR_genes_TPM_data_for_Circos_step8.1_2)
# start plotting
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_step8.1_1, plotType = "labels")
circos.initializeWithIdeogram(functional_OR_genes_TPM_data_for_Circos_step8.1_1,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.1,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(functional_OR_genes_TPM_data_for_Circos_step8.1_list, 
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      col_for_Circos_list=c("springgreen3", "lightcoral")
                      i = getI(...)
                      circos.genomicRect(region, value, 
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data.
                                         border=NA, #otherwise there are black border
                                         col = col_for_Circos_list[i])
                    })

circos.clear()


####################################################################################
   # Step 8.1 tip3: Different way 2 to plot (functional OR as example): Do not need to generate separate files for het and homo
   This is easier and better!!!!!!!! Will use this method for all the following plot.
####################################################################################
# obtain het data and make relative chr coordinates (start from 1 to sum, end = start+0.8)
functional_OR_genes_TPM_data_for_Circos_final_test <- res_functional_OR["chr_name"]
functional_OR_genes_TPM_data_for_Circos_final_test$genomic_start_relative <- NA
for (i in unique(functional_OR_genes_TPM_data_for_Circos_final_test$chr_name)) {
  functional_OR_genes_TPM_data_for_Circos_final_test[functional_OR_genes_TPM_data_for_Circos_final_test$chr_name == i, ]$genomic_start_relative <- 1:sum(functional_OR_genes_TPM_data_for_Circos_final_test$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}
functional_OR_genes_TPM_data_for_Circos_final_test$genomic_end_relative <- functional_OR_genes_TPM_data_for_Circos_final_test$genomic_start_relative + 0.8
functional_OR_genes_TPM_data_for_Circos_final_test <- cbind(functional_OR_genes_TPM_data_for_Circos_final_test, functional_OR_genes_TPM_data_for_Circos[,11:12])
functional_OR_genes_TPM_data_for_Circos_final_test <- as.data.frame(functional_OR_genes_TPM_data_for_Circos_final_test)
# start plotting
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_test, plotType = "labels")
circos.initializeWithIdeogram(functional_OR_genes_TPM_data_for_Circos_final_test,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.1,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(functional_OR_genes_TPM_data_for_Circos_final_test,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = "springgreen3")
                      circos.genomicRect(region, value,
                                         ytop.column = 2, ybottom = 0, # ytop.column=2 means the second column of value, that is the fifth of original data (homo mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = "lightcoral")
                      
                    })
circos.clear()

####################################################################################
# Step 8.1 tip4: Different way to plot functional OR with TPM above 50 set to 50
# Do not need to change TPM values (above 50) in the original file.
    This is easier and better!!!!!!!! Will use this method for all the following plot.
####################################################################################
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_final_50_as_threshold_step8.1, plotType = "labels")
circos.initializeWithIdeogram(functional_OR_genes_TPM_data_for_Circos_final_step8.1,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.1,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(functional_OR_genes_TPM_data_for_Circos_final_step8.1,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,50),
                    panel.fun = function(region, value, ...) {
                      value[[1]][value[[1]] > 50] <- 50
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border=NA, #otherwise there are black border
                                         col = ifelse(value[[2]] =="het", "springgreen3", "lightcoral")) #use [[]] instead of [] to obtain vector other than dataframe.
                    })
circos.clear()


# step 8.3: Circos plot of pseudogene OR genes --------------------------------------

#8.3 plot 1: het (green bar) and homo (red bar) use a single file with combined het and homo mean TPM
# extract the pseudogene OR data
TPM_final_for_save_step8.3 <- TPM_final_for_save_step8.1
pseudogene_OR_genes_TPM_data_for_Circos_step8.3 <- TPM_final_for_save_step8.3[grepl("^Olfr", rownames(TPM_final_for_save_step8.3)) & grepl("pseudogene", res_all$gene_type), 1:10]
pseudogene_OR_genes_TPM_data_for_Circos_step8.3$het_TPM_5samples <- rowMeans(pseudogene_OR_genes_TPM_data_for_Circos_step8.3[,1:5])
pseudogene_OR_genes_TPM_data_for_Circos_step8.3$homo_TPM_5samples <- rowMeans(pseudogene_OR_genes_TPM_data_for_Circos_step8.3[,6:10])
# obtain het data and make relative chr coordinates (start from 1 to sum, end = start+0.8)
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3 <- res_pseudogene_OR["chr_name"]
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$genomic_start_relative <- NA
for (i in unique(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$chr_name)) {
  pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3[pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$chr_name == i, ]$genomic_start_relative <- 1:sum(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$chr_name == i) #a vector from 1 to the total number of OR genes in each chromosome
}
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$genomic_end_relative <- pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$genomic_start_relative + 0.8
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3 <- cbind(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3, pseudogene_OR_genes_TPM_data_for_Circos_step8.3[,11:12])
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3 <- as.data.frame(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3)
# start plotting plot 1
pdf(paste(baseDir, "/results/pseudogene_OR_genes_TPM_Trim66_KO_Circos_plot.pdf", sep=''), paper='a4')
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3, plotType = "labels")
circos.initializeWithIdeogram(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
col_for_chr_step8.3 <- col_for_chr_universal
# obtain the colors for chr used in step8.3
col_for_chr_step8.3 <- col_for_chr_step8.3[unique(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$chr_name)]
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.3,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = "springgreen3")
                      circos.genomicRect(region, value,
                                         ytop.column = 2, ybottom = 0, # ytop.column=2 means the second column of value, that is the fifth of original data (homo mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = "lightcoral")
                      
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) \nof pseudogene OR genes in Trim66 KO",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "original TPM value", cex = 0.75)
# add legend
# check Circos mannual (https://jokergoo.github.io/circlize_book/book/index.html) and ComplexHeatmap mannual (https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html)
legend_boxs_step8.3 <- ComplexHeatmap::Legend(labels = c("het", "homo"), 
                                             type = "boxplot", 
                                             legend_gp = gpar(fill=c("springgreen3", "lightcoral")), 
                                             title_position = "topleft", 
                                             title = "Genotype")
# "npc" means Normalised Parent Coordinates (the default for unit). The origin of the viewport is (0, 0) and the viewport has a width and height of 1 unit. For example, (0.5, 0.5) is the centre of the viewport.
# check https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/unit
ComplexHeatmap::draw(legend_boxs_step8.3, x = unit(0.1, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))

circos.clear()

# 8.3 plot 2: since 98%, and 99% percentil of the TPM values is 4.79 and 6.07, I will manipulate the TPM value above 10 to be 10 for better visualizaiton.
quantile(c(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$het_TPM_5samples, 
           pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$homo_TPM_5samples),
         c(0.5, 0.75, 0.9,0.95, 0.98,0.99,0.995,1))
sum(c(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$het_TPM_5samples, 
      pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3$homo_TPM_5samples) >10)  #total 4 data points with TPM above 10
pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.3 <- pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3
pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.3[pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.3$het_TPM_5samples>10, ]$het_TPM_5samples <- 10
pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.3[pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.3$homo_TPM_5samples>10, ]$homo_TPM_5samples <- 10
# start plotting plot 2
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_step8.1_final_40_as_threshold, plotType = "labels")
circos.initializeWithIdeogram(pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.3,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.3,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.3,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = "springgreen3")
                      circos.genomicRect(region, value,
                                         ytop.column = 2, ybottom = 0, # ytop.column=2 means the second column of value, that is the fifth of original data (homo mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = "lightcoral")
                      
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) \nof pseudogene OR genes in Trim66 KO",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "TPM larger than 10 are \nchanged to 10 for better visualization\n(4 data points out of 271*2 are changed)", cex = 0.75)
# add legend
ComplexHeatmap::draw(legend_boxs_step8.3, x = unit(0.1, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
circos.clear()
dev.off()

# step 8.4: Circos plot of pseudogene OR genes with significant changes marked--------------------------------------

#8.4 plot 1: het (green bar) and homo (red bar) use a single file with combined het and homo mean TPM
# extract the pseudogene OR data
pseudogene_OR_genes_TPM_data_for_Circos_step8.4 <- pseudogene_OR_genes_TPM_data_for_Circos_step8.3
# obtain het data and make relative chr coordinates (start from 1 to sum, end = start+0.8)
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4 <- pseudogene_OR_genes_TPM_data_for_Circos_final_step8.3
# add one column to indicate if significantly changed
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4$is_significant_padj_0.05_log2FC_0.585 <- res_pseudogene_OR$padj < 0.05 & abs(res_pseudogene_OR$log2FoldChange)>=0.585
# change NA to FALSE
pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4[is.na(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4$is_significant_padj_0.05_log2FC_0.585), ]$is_significant_padj_0.05_log2FC_0.585 <- FALSE
# start plotting plot 1
pdf(paste(baseDir, "/results/pseudogene_OR_genes_TPM_Trim66_KO_Circos_plot_sig_nonsig_labeled.pdf", sep=''), paper='a4')
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(pseudogene_OR_genes_TPM_data_for_Circos_step8.4_step8.1_final_40_as_threshold, plotType = "labels")
circos.initializeWithIdeogram(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
col_for_chr_step8.4 <- col_for_chr_step8.3
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.4,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = ifelse(value[[3]] ==TRUE, "springgreen3", "gray80"))
                      circos.genomicRect(region, value,
                                         ytop.column = 2, ybottom = 0, # ytop.column=2 means the second column of value, that is the fifth of original data (homo mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = ifelse(value[[3]] ==TRUE, "lightcoral", "gray10"))
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) \nof pseudogene OR genes in Trim66 KO",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "original TPM value", cex = 0.75)
# add legend
# check Circos mannual (https://jokergoo.github.io/circlize_book/book/index.html) and ComplexHeatmap mannual (https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html)
legend_boxs_step8.4 <- ComplexHeatmap::Legend(labels = c("het sig", "homo sig", "het non_sig", "homo_non_sig"), 
                                              type = "boxplot", 
                                              legend_gp = gpar(fill=c("springgreen3", "lightcoral", "gray80", "gray10")), 
                                              title_position = "topleft", 
                                              title = "Genotype")
# "npc" means Normalised Parent Coordinates (the default for unit). The origin of the viewport is (0, 0) and the viewport has a width and height of 1 unit. For example, (0.5, 0.5) is the centre of the viewport.
# check https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/unit
ComplexHeatmap::draw(legend_boxs_step8.4, x = unit(0.1, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))

circos.clear()

# 8.4 plot 2: since 98%, and 99% percentil of the TPM values is 4.79 and 6.07, I will manipulate the TPM value above 10 to be 10 for better visualizaiton.
quantile(c(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4$het_TPM_5samples, 
           pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4$homo_TPM_5samples),
         c(0.5, 0.75, 0.9,0.95, 0.98,0.99,0.995,1))
sum(c(pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4$het_TPM_5samples, 
      pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4$homo_TPM_5samples) >10)  #total 4 data points with TPM above 10
pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.4 <- pseudogene_OR_genes_TPM_data_for_Circos_final_step8.4
pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.4[pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.4$het_TPM_5samples>10, ]$het_TPM_5samples <- 10
pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.4[pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.4$homo_TPM_5samples>10, ]$homo_TPM_5samples <- 10
# start plotting plot 2
circos.par(start.degree=90, 
           gap.after=0.1, # gap.after to adjust the gap between different sectors
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(functional_OR_genes_TPM_data_for_Circos_step8.1_final_40_as_threshold, plotType = "labels")
circos.initializeWithIdeogram(pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.4,
                              plotType = "labels")  # plotType = c("axis","labels") to show axis
# add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step8.4,
             bg.lwd =0.1, # adjus the thickness of border
             bg.border = "gray80", track.height = 0.05) #track.height is the relative height to the radius of 1
# add the second track of real data
circos.genomicTrack(pseudogene_OR_genes_TPM_data_for_Circos_final_10_as_threshold_step8.4,
                    bg.border = NA, # otherwise there are borders of different sectors
                    track.height = 0.5, # set the track height, which is is the percentage to the radius of the unit circles
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop.column = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = ifelse(value[[3]] ==TRUE, "springgreen3", "gray80"))
                      circos.genomicRect(region, value,
                                         ytop.column = 2, ybottom = 0, # ytop.column=2 means the second column of value, that is the fifth of original data (homo mean TPM).
                                         border=NA, #otherwise there are black border
                                         col = ifelse(value[[3]] ==TRUE, "lightcoral", "gray10"))
                    })
title(main = "Circos plot for TPM (mean of 5 het or homo samples) of pseudogene OR genes in Trim66 KO \nwith information of signicant changes padj<0.05 and abs(log2FC)>=0.585",
      font.main=2, cex.main=0.8, line = -1) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "TPM larger than 10 are \nchanged to 10 for better visualization\n(4 data points out of 271*2 are changed)", cex = 0.75)
# add legend
ComplexHeatmap::draw(legend_boxs_step8.4, x = unit(0.1, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
circos.clear()
dev.off()


# Step 9: GO_analysis -------------------------------------------------------------
#####################################################################################################################
          #Plot step 9: GO analysis
#####################################################################################################################
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(pathview)
library(GOSemSim)
library(org.Mm.eg.db)
library(rJava)
library(RDAVIDWebService)

res_homo_vs_het_sig_step9 <- res_all[which(res_all$padj<0.05 & abs(res_all$log2FoldChange)>=0.585),]
res_homo_vs_het_sig_step9@rownames
# bitr to convert the gene symbol to other types that can be used for GO analysis
res_homo_vs_het_sig_genes_step9 <- bitr(res_homo_vs_het_sig_step9@rownames, fromType = "SYMBOL",
                                        toType = c("ENSEMBL", "ENTREZID"),
                                        OrgDb = org.Mm.eg.db)

# set up the background gene list of MOE specific genes with the criteria of baseMean>=20
MOE_specific_genes_step9 <- res_all[which(res_all$baseMean >= 20), ]
background_MOE_specific_genes_step9 <- bitr(MOE_specific_genes_step9@rownames, fromType = "SYMBOL",
                                            toType = c("ENSEMBL", "ENTREZID"),
                                            OrgDb = org.Mm.eg.db) 
# 9.1 Reactome enrichment
reactome_pathway_step9 <- enrichPathway(res_homo_vs_het_sig_genes_step9$ENTREZID, organism = "mouse")
head(reactome_pathway_step9, 5)
clusterProfiler::dotplot(reactome_pathway_step9)
barplot(reactome_pathway_step9, drop=TRUE, showCategory=12)
# 9.2 GO enrichment ont = One of "MF", "BP", and "CC" subontologies. MF: molecular functions; BP: biological pathway; CC: cellular component
ggo_homo_vs_het_step9 <- groupGO(gene     = res_homo_vs_het_sig_genes_step9$ENTREZID,
                                 OrgDb    = org.Mm.eg.db,
                                 ont      = "CC",
                                 level    = 5,  #https://www.biostars.org/p/237816/#237840 for levels
                                 readable = TRUE)
head(ggo_homo_vs_het_step9)
barplot(ggo_homo_vs_het_step9, drop=TRUE, showCategory=12)
clusterProfiler::dotplot(ggo_homo_vs_het_step9, showCategory=12)

# GO over-representation test. In my understanding, GO over-representation and enrichment test are similar thing, just using different statistical methods
ego_homo_vs_het_step9 <- clusterProfiler::enrichGO(gene          = res_homo_vs_het_sig_genes_step9$ENTREZID,
                                                   OrgDb         = org.Mm.eg.db,
                                                   ont           = "MF",
                                                   pAdjustMethod = "BH",
                                                   universe = background_MOE_specific_genes_step9$ENTREZID,
                                                   minGSSize = 10,
                                                   maxGSSize = 500,
                                                   pvalueCutoff  = 0.05,
                                                   qvalueCutoff  = 0.05,
                                                   readable      = TRUE)
head(ego_homo_vs_het_step9, 20)
barplot(ego_homo_vs_het_step9, drop=TRUE, showCategory=20)
pdf(paste(baseDir, "/results/GO_analysis_Trim66_KO_padj_0.05_log2FC_0.585.pdf", sep=''), 
    paper = 'a4')
dotplot(ego_homo_vs_het_step9, showCategory=20) #used this as GO enrichment for our data
dev.off()

# 9.3 KEGG enrichment
# KEGG over-representation test
kk_homo_vs_het_step9 <- enrichKEGG(gene         = res_homo_vs_het_sig_genes_step9$ENTREZID,
                                   organism     = 'mmu',
                                   pvalueCutoff = 0.05)
head(kk_homo_vs_het_step9)
dotplot(kk_homo_vs_het_step9, showCategory=20)
barplot(kk_homo_vs_het_step9, drop=TRUE, showCategory=20)

# KEGG Module over-representation test
mkk_homo_vs_het_step9 <- enrichMKEGG(gene = res_homo_vs_het_sig_genes_step9$ENTREZID,
                                     organism = 'mmu')
head(mkk_homo_vs_het_step9)
dotplot(mkk_homo_vs_het_step9)
barplot(mkk_homo_vs_het_step9, drop=TRUE, showCategory=12)

# 9.4 DAVID functional analysis
david_homo_vs_het_step9 <- enrichDAVID(gene = res_homo_vs_het_sig_genes_step9$ENTREZID,
                                       idType = "ENTREZ_GENE_ID",
                                       annotation = "GOTERM_BP_FAT",
                                       david.user = "clusterProfiler@hku.hk")
head(david_homo_vs_het_step9)
dotplot(david_homo_vs_het_step9)
barplot(david_homo_vs_het_step9, drop=TRUE, showCategory=12)


# Step 10: Circos plot of all genes and show OR clusters ------------------

#####################################################################################################################
            # Plot Step 10: plot Circos plot for all genes and show OR clusters
#####################################################################################################################
#plot step10: 
# circos plot with track 1: chr with different colors;
# track 2: all of the upregulated (lightcoral color) and downregulated (springgreen3 color) and unchanged (gray80 color) genes in homo. Also add the information of functional 
# track 3: OR clusters and TAAR cluster, indicate TAAR cluster with blue color;
# track 4: position of all functional OR and TAAR genes, indicate TAARs with blue color;
# track 5: position of all pseudogene OR genes.
#########################################################################################
################################### prepare the data for circos plotting ################
res_all_step10 <- res_all_step3.5 # has the column to indicate if significant changes, threshold: padj<0.05 and abs(log2FC)>=0.585
res_all_dataframe_step10 <- as.data.frame(res_all_step10)
res_all_dataframe_for_Circos_step10 <-res_all_dataframe_step10[, c("chr_name", "genomic_start", "genomic_end", "log2FoldChange", "if_significant_upregulate_or_downregulate")]
# change NA in log2FoldChange to 0
sum(is.na(res_all_dataframe_for_Circos_step10$log2FoldChange))
res_all_dataframe_for_Circos_step10[is.na(res_all_dataframe_for_Circos_step10$log2FoldChange), ]$log2FoldChange <- 0
sum(is.na(res_all_dataframe_for_Circos_step10$log2FoldChange))
sum(res_all_dataframe_for_Circos_step10$log2FoldChange == 0)

################################################################
                   # obtain OR cluster information
################################################################
# make a varible to store OR cluster information
res_all_ORs_step10 <- res_all[grepl("^Olfr", rownames(res_all)), ]
dim(res_all_ORs_step10)
# add one column to extract chromosome numbers, chrX to 20. There are no ORs in chr18 and chr20.
unique(res_all_ORs_step10$chr_name) # check the order of chr in unique(res_all_ORs_step10$chr_name). You can use gtools::mixedsort(unique(res_all_ORs_step10$chr_name)) and make "chr_name" column as factors if it is ordered as chr1, chr10, chr11,....
res_all_ORs_step10$chr_number_only <- ifelse(res_all_ORs_step10$chr_name == "chrX", 20, as.numeric(substring(res_all_ORs_step10$chr_name, 4))) #check substring function, https://stackoverflow.com/questions/17215789/extract-a-substring-according-to-a-pattern
# add one column to show if the two ORs are above 1 Mb apart. Longzhi said :因为我画过or间距（log10之后）的histogram，1mb那里是个谷
res_all_ORs_step10$above_1Mb_between_consective_ORs <- diff(c(1, res_all_ORs_step10$genomic_start)) >= 1000000 &
  diff(c(1, res_all_ORs_step10$chr_number_only)) == 0 # add a row with number 1 for diff function. Also make sure it is in the same chromosome
# change the "above_1Mb_between_consective_ORs" value of the first chr as TRUE
for (i in unique(res_all_ORs_step10$chr_name)) {
  res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$above_1Mb_between_consective_ORs[1] <- TRUE #the first as TRUE
}
# add two columns to indicate each and total OR cluster number in each chromosome with initial value as 0 
res_all_ORs_step10$total_number_OR_clusters_in_each_chr_with_chr_name <- ""
res_all_ORs_step10$total_number_OR_clusters_in_each_chr <- 0
res_all_ORs_step10$number_each_OR_cluster_in_each_chr_with_chr_name <- ""
res_all_ORs_step10$number_each_OR_cluster_in_each_chr <- 0
# change the above two columns with real data.
for (i in unique(res_all_ORs_step10$chr_name)) {
  res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$total_number_OR_clusters_in_each_chr_with_chr_name <- paste("total_", sum(res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$above_1Mb_between_consective_ORs), "_OR_clusters_in_", i, sep = "")
  res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$total_number_OR_clusters_in_each_chr <- sum(res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$above_1Mb_between_consective_ORs)
  res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$number_each_OR_cluster_in_each_chr_with_chr_name <- paste("OR_cluster_No.", cumsum(res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$above_1Mb_between_consective_ORs), "_in_", i, sep = "")
  res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$number_each_OR_cluster_in_each_chr <- cumsum(res_all_ORs_step10[res_all_ORs_step10$chr_name == i, ]$above_1Mb_between_consective_ORs)
}
# add two columns to indicate each and total OR cluster number in the whole chromosome. 
res_all_ORs_step10$total_number_OR_clusters_in_whole_chr_with_name <- paste("total_", sum(res_all_ORs_step10$above_1Mb_between_consective_ORs), "_OR_clusters_in_whole_genome", sep = "")
res_all_ORs_step10$total_number_OR_clusters_in_whole_chr <- sum(res_all_ORs_step10$above_1Mb_between_consective_ORs)
res_all_ORs_step10$number_each_OR_cluster_in_whole_chr_with_chr_name <- paste("OR_cluster_No.", cumsum(res_all_ORs_step10$above_1Mb_between_consective_ORs), "_of_whole_genome_in_", res_all_ORs_step10$chr_name, sep = "")
res_all_ORs_step10$number_each_OR_cluster_in_whole_chr <- cumsum(res_all_ORs_step10$above_1Mb_between_consective_ORs)
head(res_all_ORs_step10)
colnames(res_all_ORs_step10)
# convert to dataframe
res_all_ORs_dataframe_step10 <- as.data.frame(res_all_ORs_step10)
# save this information of all ORs with OR cluster IDs to file in common database
write.csv(res_all_ORs_dataframe_step10[, -(1:6)], file = "/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/All_ORs_with_OR_clusters_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25.csv")

# create a dataframe to store chr_name, genomic start, genomic end, and other information of each OR cluster
OR_clusters_data_for_Circos_step10 <- data.frame(matrix(nrow = sum(res_all_ORs_dataframe_step10$above_1Mb_between_consective_ORs), ncol = 13))
colnames(OR_clusters_data_for_Circos_step10) <- c("chr_name",
                                                  "genomic_start",
                                                  "genomic_end",
                                                  "OR_cluster_length",
                                                  "total_number_OR_clusters_in_whole_chr_with_name",
                                                  "total_number_OR_clusters_in_whole_chr",
                                                  "total_number_OR_clusters_in_each_chr_with_chr_name",
                                                  "total_number_OR_clusters_in_each_chr",
                                                  "number_this_OR_cluster_in_each_chr_with_chr_name",
                                                  "number_this_OR_cluster_in_each_chr",
                                                  "how_many_functional_ORs_in_this_cluster",
                                                  "how_many_pseudogene_ORs_in_this_cluster",
                                                  "how_many_all_ORs_in_this_cluster")
rownames(OR_clusters_data_for_Circos_step10) <- unique(res_all_ORs_dataframe_step10$number_each_OR_cluster_in_whole_chr_with_chr_name)
for (i in rownames(OR_clusters_data_for_Circos_step10)) {
  res_all_ORs_temp_step10 <- res_all_ORs_dataframe_step10[res_all_ORs_dataframe_step10$number_each_OR_cluster_in_whole_chr_with_chr_name == i, ]
  OR_clusters_data_for_Circos_step10[i, c(1:3, 5:10)] <- res_all_ORs_temp_step10[1,c("chr_name",
                                                                                     "genomic_start",
                                                                                     "genomic_end",
                                                                                     "total_number_OR_clusters_in_whole_chr_with_name",
                                                                                     "total_number_OR_clusters_in_whole_chr",
                                                                                     "total_number_OR_clusters_in_each_chr_with_chr_name",
                                                                                     "total_number_OR_clusters_in_each_chr",
                                                                                     "number_each_OR_cluster_in_each_chr_with_chr_name",
                                                                                     "number_each_OR_cluster_in_each_chr")]
  OR_clusters_data_for_Circos_step10[i, 3] <- res_all_ORs_temp_step10[nrow(res_all_ORs_temp_step10), "genomic_end"] #change the genomic_end to the last element of res_all_ORs_temp_step10
  OR_clusters_data_for_Circos_step10[i, 4] <- OR_clusters_data_for_Circos_step10[i,3]-OR_clusters_data_for_Circos_step10[i,2]
  OR_clusters_data_for_Circos_step10[i, 11] <- sum(grepl("protein_coding", res_all_ORs_temp_step10$gene_type))
  OR_clusters_data_for_Circos_step10[i, 12] <- sum(grepl("pseudogene", res_all_ORs_temp_step10$gene_type))
  OR_clusters_data_for_Circos_step10[i, 13] <- nrow(res_all_ORs_temp_step10)
}
# check if the information is correct
head(OR_clusters_data_for_Circos_step10)
tail(OR_clusters_data_for_Circos_step10)
sum(OR_clusters_data_for_Circos_step10$how_many_functional_ORs_in_this_cluster)
sum(OR_clusters_data_for_Circos_step10$how_many_pseudogene_ORs_in_this_cluster)
sum(OR_clusters_data_for_Circos_step10$how_many_all_ORs_in_this_cluster)
sum(OR_clusters_data_for_Circos_step10$how_many_functional_ORs_in_this_cluster) == nrow(res_functional_OR)
sum(OR_clusters_data_for_Circos_step10$how_many_pseudogene_ORs_in_this_cluster) == nrow(res_pseudogene_OR)
# add two columns to show how many large OR clusters containing >= 10 functional ORs in whole chromosomes
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10functional_ORs_in_whole_chr_with_name <- paste("total_", 
                                                                                                                             sum(OR_clusters_data_for_Circos_step10$how_many_functional_ORs_in_this_cluster >=10), 
                                                                                                                             "_large_OR_clusters(>=10_functional_ORs)_in_whole_genome", 
                                                                                                                             sep = "")
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10functional_ORs_in_whole_chr <- sum(OR_clusters_data_for_Circos_step10$how_many_functional_ORs_in_this_cluster >=10)
# add two columns to show how many large OR clusters containing >= 10 functional ORs in each chromosome
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10functional_ORs_in_each_chr_with_chr_name <- ""
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10functional_ORs_in_each_chr <- 0
for (i in OR_clusters_data_for_Circos_step10$chr_name) {
  OR_clusters_data_for_Circos_temp_step10 <- OR_clusters_data_for_Circos_step10[OR_clusters_data_for_Circos_step10$chr_name == i, ]
  OR_clusters_data_for_Circos_step10[OR_clusters_data_for_Circos_step10$chr_name == i, ]$total_number_large_OR_clusters_more_than_10functional_ORs_in_each_chr_with_chr_name <- paste("total_", 
                                                                                                                                                                                      sum(OR_clusters_data_for_Circos_temp_step10$how_many_functional_ORs_in_this_cluster >=10), 
                                                                                                                                                                                      "_large_OR_clusters(>=10_functional_ORs)_in_", 
                                                                                                                                                                                      i, 
                                                                                                                                                                                      sep = "")
  OR_clusters_data_for_Circos_step10[OR_clusters_data_for_Circos_step10$chr_name == i, ]$total_number_large_OR_clusters_more_than_10functional_ORs_in_each_chr <- sum(OR_clusters_data_for_Circos_temp_step10$how_many_functional_ORs_in_this_cluster >=10)
}



# add two columns to show how many large OR clusters containing >= 10 all ORs (including functional and pseudogene) in whole chromosomes
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10all_ORs_in_whole_chr_with_name <- paste("total_", 
                                                                                                                      sum(OR_clusters_data_for_Circos_step10$how_many_all_ORs_in_this_cluster >=10), 
                                                                                                                      "_large_OR_clusters(>=10_functional_ORs)_in_whole_genome", 
                                                                                                                      sep = "")
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10all_ORs_in_whole_chr <- sum(OR_clusters_data_for_Circos_step10$how_many_all_ORs_in_this_cluster >=10)
# and two more columns to show how many large OR clusters containing >= 10 all ORs (including functional and pseudogene) in each chromosome
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10all_ORs_in_each_chr_with_chr_name <- ""
OR_clusters_data_for_Circos_step10$total_number_large_OR_clusters_more_than_10all_ORs_in_each_chr <- 0
for (i in OR_clusters_data_for_Circos_step10$chr_name) {
  OR_clusters_data_for_Circos_temp_step10 <- OR_clusters_data_for_Circos_step10[OR_clusters_data_for_Circos_step10$chr_name == i, ]
  OR_clusters_data_for_Circos_step10[OR_clusters_data_for_Circos_step10$chr_name == i, ]$total_number_large_OR_clusters_more_than_10all_ORs_in_each_chr_with_chr_name <- paste("total_", 
                                                                                                                                                                               sum(OR_clusters_data_for_Circos_temp_step10$how_many_all_ORs_in_this_cluster >=10), 
                                                                                                                                                                               "_large_OR_clusters(>=10_all_ORs)_in_", 
                                                                                                                                                                               i, 
                                                                                                                                                                               sep = "")
  OR_clusters_data_for_Circos_step10[OR_clusters_data_for_Circos_step10$chr_name == i, ]$total_number_large_OR_clusters_more_than_10all_ORs_in_each_chr <- sum(OR_clusters_data_for_Circos_temp_step10$how_many_all_ORs_in_this_cluster >=10)
}

# save to file in common database
write.csv(OR_clusters_data_for_Circos_step10, file = "/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/OR_clusters_full_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25.csv")

# extract data from OR_clusters_data_for_Circos_step10 to obtain tidy version with summary of OR clusters
OR_clusters_data_summary_step10 <- data.frame(matrix(nrow = length(unique(OR_clusters_data_for_Circos_step10$chr_name)), ncol = 19))
colnames(OR_clusters_data_summary_step10) <- colnames(OR_clusters_data_for_Circos_step10)[-(9:10)]
rownames(OR_clusters_data_summary_step10) <- unique(res_all_ORs_dataframe_step10$chr_name)
for (i in rownames(OR_clusters_data_summary_step10)) {
  res_all_ORs_temp_step10 <- OR_clusters_data_for_Circos_step10[OR_clusters_data_for_Circos_step10$chr_name == i, ]
  OR_clusters_data_summary_step10[i, ] <- res_all_ORs_temp_step10[1, colnames(OR_clusters_data_summary_step10)]
  OR_clusters_data_summary_step10[i, "genomic_end"] <- res_all_ORs_temp_step10[nrow(res_all_ORs_temp_step10), "genomic_end"]
  OR_clusters_data_summary_step10[i, "OR_cluster_length"] <- OR_clusters_data_summary_step10[i, "genomic_end"] - OR_clusters_data_summary_step10[i, "genomic_start"]
  OR_clusters_data_summary_step10[i, "how_many_functional_ORs_in_this_cluster"] <- sum(res_all_ORs_temp_step10$how_many_functional_ORs_in_this_cluster)
  OR_clusters_data_summary_step10[i, "how_many_pseudogene_ORs_in_this_cluster"] <- sum(res_all_ORs_temp_step10$how_many_pseudogene_ORs_in_this_cluster)
  OR_clusters_data_summary_step10[i, "how_many_all_ORs_in_this_cluster"] <- sum(res_all_ORs_temp_step10$how_many_all_ORs_in_this_cluster)
}
colnames(OR_clusters_data_summary_step10)[4] <- "ORs_in_each_chr_range"
# check if the information is correct
OR_clusters_data_summary_step10
sum(OR_clusters_data_summary_step10$total_number_OR_clusters_in_each_chr) == unique(OR_clusters_data_summary_step10$total_number_OR_clusters_in_whole_chr)
sum(OR_clusters_data_summary_step10$total_number_large_OR_clusters_more_than_10functional_ORs_in_each_chr) == unique(OR_clusters_data_summary_step10$total_number_large_OR_clusters_more_than_10functional_ORs_in_whole_chr)
sum(OR_clusters_data_summary_step10$total_number_large_OR_clusters_more_than_10all_ORs_in_each_chr) == unique(OR_clusters_data_summary_step10$total_number_large_OR_clusters_more_than_10all_ORs_in_whole_chr)
# save to file in common database
write.csv(OR_clusters_data_summary_step10, file = "/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/OR_clusters_summary_info_by_my_criteria_OR_distance_above_1Mb_Gencode_vM25.csv")

################################################################
              # obtain TAAR cluster information
################################################################
# make a separate file for TAAR information
res_all_TAARs_step10 <- res_all[grepl("^Taar", rownames(res_all)), ]
dim(res_all_TAARs_step10)
TAAR_cluster_data_for_Circos_step10 <- data.frame(matrix(nrow = 1, ncol = 8))
colnames(TAAR_cluster_data_for_Circos_step10) <- c("chr_name", 
                                                   "genomic_start",
                                                   "genomic_end",
                                                   "TAAR_cluster_length",
                                                   "how_many_total_TAARs",
                                                   "how_many_functional_TAARs_including_TAAR1",
                                                   "how_many_functional_TAARs_excluding_TAAR1",
                                                   "how_many_pseudogene_TAARs")

TAAR_cluster_data_for_Circos_step10[1,"chr_name"] <- res_all_TAARs_step10["Taar1", "chr_name"]
TAAR_cluster_data_for_Circos_step10[1,"genomic_start"] <- res_all_TAARs_step10["Taar1", "genomic_start"]
TAAR_cluster_data_for_Circos_step10[1,"genomic_end"] <- res_all_TAARs_step10["Taar9", "genomic_start"]
TAAR_cluster_data_for_Circos_step10[1,"TAAR_cluster_length"] <- TAAR_cluster_data_for_Circos_step10[1,"genomic_end"] - TAAR_cluster_data_for_Circos_step10[1,"genomic_start"]
TAAR_cluster_data_for_Circos_step10[1,"how_many_total_TAARs"] <- nrow(res_all_TAARs_step10)
TAAR_cluster_data_for_Circos_step10[1,"how_many_functional_TAARs_including_TAAR1"] <- sum(grepl("protein_coding", res_all_TAARs_step10$gene_type))
TAAR_cluster_data_for_Circos_step10[1,"how_many_functional_TAARs_excluding_TAAR1"] <- sum(grepl("^Taar[2-9]", rownames(res_all_TAARs_step10)) & grepl("protein_coding", res_all_TAARs_step10$gene_type))
TAAR_cluster_data_for_Circos_step10[1,"how_many_pseudogene_TAARs"] <- sum(grepl("pseudogene", res_all_TAARs_step10$gene_type))
# save to file in common database
write.csv(TAAR_cluster_data_for_Circos_step10, file = "/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/TAAR_cluster_summary_info_Gencode_vM25.csv")


#########################################################################################
################################### Start Circos plotting #############################
# start plotting
pdf(paste(baseDir, "/results/all_genes_Trim66_KO_Circos_plot.pdf", sep=''), paper='a4')
# change the basic settings of sectors and tracks
circos.par(start.degree=90, 
           gap.after=0, # gap.after to adjust the gap between different sectors
           track.margin=c(0.008,0.008), # set the bottom and topmargin between tracks
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(res_all_dataframe_for_Circos_step10, plotType = "labels")
circos.initializeWithIdeogram(res_all_dataframe_for_Circos_step10,
                              plotType = c("labels"))  # plotType = c("axis","labels") to show axis and labels
# 1. add the first track of the different chromosomes with different colors
col_for_chr_step10 <- col_for_chr_universal
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step10,
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the second track of log2FoldChange of all genes with position of OR and TAAR clusters as background bars
OR_and_TAAR_clusters_and_all_gene_changes_for_Circos_list_step10 = list(OR_clusters_data_for_Circos_step10,
                                                                        TAAR_cluster_data_for_Circos_step10,
                                                                        res_all_dataframe_for_Circos_step10)
circos.genomicTrack(OR_and_TAAR_clusters_and_all_gene_changes_for_Circos_list_step10,
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2, # line type, check https://r-charts.com/base-r/line-types/
                    bg.border="gray95",
                    track.height = 0.4, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(-8, 10), # range of -7.98 to 9.75 by range(res_all_dataframe_for_Circos_step10$log2FoldChange)
                    panel.fun = function(region, value, ...) {
                      col_for_Circos_list_step10=c("gold4", "blue")
                      tranparency_list_step10=c(0.4, 0) #vector to adjust the transparency of color
                      i = getI(...) # i is the index for the current data frame from the list. getI(...) has to be used inside the panel.fun.
                      if (i < 3) {
                        circos.genomicRect(region, value, # plot OR/TAAR clusters first as background
                                           ytop = 10, ybottom = -8, # ybottom =-18 so that the rectangele can extend to next track (ylim is -8 to 10 in the next track).
                                           border = add_transparency(col_for_Circos_list_step10[i], 
                                                                     transparency = tranparency_list_step10[i]), #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency.
                                           col = add_transparency(col_for_Circos_list_step10[i], 
                                                                  transparency = tranparency_list_step10[i])) #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency.
                      } else {
                        circos.genomicPoints(region, value, # plot genes with no signicant changes firstly
                                             numeric.column = 1,  # numeric.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                             cex = 0.4, pch = 16,
                                             col = add_transparency(hutils::Switch(value[[2]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                                   Not_changed="gray80",
                                                                                   Downregulated_significantly="springgreen3", 
                                                                                   Upregulated_significantly="lightcoral",
                                                                                   DEFAULT = "white"), 
                                                                    transparency = ifelse(value[[2]] == "Not_changed", 0, 1))) #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency. check https://cran.r-project.org/web/packages/circlize/circlize.pdf
                        circos.genomicPoints(region, value,
                                             numeric.column = 1,  # numeric.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                             cex = 0.4, pch = 16,
                                             col = add_transparency(hutils::Switch(value[[2]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                                   Downregulated_significantly="springgreen3", 
                                                                                   Upregulated_significantly="lightcoral",
                                                                                   DEFAULT = "white"), 
                                                                    transparency = ifelse(value[[2]] == "Not_changed", 1, 0))) #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency. check https://cran.r-project.org/web/packages/circlize/circlize.pdf
                        circos.lines(CELL_META$cell.xlim, c(0, 0), # plot the line at y=0
                                     lty = 1, 
                                     lwd = 0.25,
                                     col = "gray50")
                      }
                    })
# add text for legend indicating this track
circos.text(CELL_META$xcenter, 
            (CELL_META$cell.ylim[2])/2,
            "Track 1: Trim66 homo vs het with OR and TAAR clusters shown",
            sector.index = "chrY",
            facing = "bending.inside",
            cex = 0.25)

# 3. add the third track of position of functional OR genes (gray) and TAAR genes (blue)
res_functional_OR_and_TAAR_step10 <- res_functional_OR_and_TAAR_step3.3 #has the column above_1Mb_between_consective_ORs
res_functional_OR_and_TAAR_dataframe_step10 <- as.data.frame(res_functional_OR_and_TAAR_step10)
res_functional_OR_and_TAAR_dataframe_for_Circos_step10 <- res_functional_OR_and_TAAR_dataframe_step10[, c("chr_name", "genomic_start", "genomic_end")]
# add one column to indicate TAAR or OR
res_functional_OR_and_TAAR_dataframe_for_Circos_step10$is_TAAR_or_OR <- "OR"
res_functional_OR_and_TAAR_dataframe_for_Circos_step10[grepl("Taar", rownames(res_functional_OR_and_TAAR_dataframe_for_Circos_step10)), ]$is_TAAR_or_OR <- "TAAR"
circos.genomicTrack(res_functional_OR_and_TAAR_dataframe_for_Circos_step10,
                    bg.lwd = 0.01, # adjust the thickness of border between sectors
                    bg.lty = 2,
                    bg.border="gray95",
                    track.height = 0.1, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border=ifelse(value[[1]] == "OR", "gray10", "blue"),
                                         col = ifelse(value[[1]] == "OR", "gray10", "blue")) #use [[]] instead of [] to obtain vector other than dataframe.
                    })
# add text for legend indicating this track
circos.text(CELL_META$xcenter, 
            CELL_META$ycenter,
            "Track 2: Functional ORs and TAARs",
            sector.index = "chrY",
            facing = "bending.inside",
            cex = 0.25)
# 4. add the fourth track of position of pseudogene OR genes
res_pseudogene_OR_step10 <- res_pseudogene_OR_step3.2 #has the column above_1Mb_between_consective_ORs
res_pseudogene_OR_dataframe_step10 <- as.data.frame(res_pseudogene_OR_step10)
res_pseudogene_OR_dataframe_for_Circos_step10 <- res_pseudogene_OR_dataframe_step10[, c("chr_name", "genomic_start", "genomic_end")]
circos.genomicTrack(res_pseudogene_OR_dataframe_for_Circos_step10,
                    bg.lwd = 0.01, # otherwise there are borders of different sectors
                    bg.lty = 2,
                    bg.border="gray95",
                    track.height = 0.1, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop = 1, ybottom = 0, # ytop.column=1 means the first column of value, that is the fourth of original functional_OR_genes_TPM_data_for_Circos_final data.
                                         border="gray30", #otherwise there are black border
                                         col = "gray30") #use [[]] instead of [] to obtain vector other than dataframe.
                    })
# add text for legend indicating this track
circos.text(CELL_META$xcenter, 
            CELL_META$ycenter,
            "Track 3: OR pseudogenes",
            sector.index = "chrY",
            facing = "bending.inside",
            cex = 0.25)
# add title
title(main = "Circos plot for changes of all genes (homo vs het) in Trim66 KO",
      font.main=2, cex.main=0.8, line = -0.5) #cex.main means the magnification for main title. line means the distance between title and plot border
text(0, 0, "Significant: padj<0.05 & \nabs(log2FC)>=0.585 \nblue bar: TAAR cluster", cex = 0.6)
# add legend
# check Circos mannual (https://jokergoo.github.io/circlize_book/book/index.html) and ComplexHeatmap mannual (https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html)
legend_boxs_step10 <- ComplexHeatmap::Legend(labels = c("Significant downreuglation", "Significant upregulation", "no sigificant change"), 
                                             type = "points", 
                                             legend_gp = gpar(col=c("springgreen3", "lightcoral", "gray80")), 
                                             title_position = "topleft", 
                                             title = "Genotype")
# "npc" means Normalised Parent Coordinates (the default for unit). The origin of the viewport is (0, 0) and the viewport has a width and height of 1 unit. For example, (0.5, 0.5) is the centre of the viewport.
# check https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/unit
ComplexHeatmap::draw(legend_boxs_step10, x = unit(0.05, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
circos.clear()
dev.off()


##################################################################################
      # Summary on plot step 10: To plot OR/TAAR clusters and gene changes separately
##################################################################################
# Summary 1. add positions of OR and TAAR clusters.
OR_and_TAAR_clusters_for_Circos_list_step10 = list(OR_clusters_data_for_Circos_step10,
                                                   TAAR_cluster_data_for_Circos_step10)
circos.genomicTrack(OR_and_TAAR_clusters_for_Circos_list_step10, 
                    bg.lwd = 0.2, # adjust the thickness of border between sectors
                    bg.lty = "1F",
                    bg.border="gray80",
                    track.height = 0.1, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,1), 
                    panel.fun = function(region, value, ...) {
                      col_for_Circos_list_step10=c("black", "blue")
                      tranparency_list_step10=c(0.1, 0)
                      i = getI(...)
                      circos.genomicRect(region, value, 
                                         ytop = 1, ybottom = 0, # ybottom =-18 so that the rectangele can extend to next track (ylim is -8 to 10 in the next track).
                                         border= add_transparency(col_for_Circos_list_step10[i], 
                                                                  transparency = tranparency_list_step10[i]), #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency.
                                         col = add_transparency(col_for_Circos_list_step10[i], 
                                                                transparency = tranparency_list_step10[i])) #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency.
                    })
# Summary 2. add the second track of log2FoldChange of all genes with position of OR and TAAR clusters as background bars
circos.genomicTrack(res_all_dataframe_for_Circos_step10,
                    bg.lwd = 0.2, # adjust the thickness of border between sectors
                    bg.lty = "1F",
                    bg.border="gray80",
                    track.height = 0.3, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(-8, 10), # range of -7.98 to 9.75 by range(res_all_dataframe_for_Circos_step10$log2FoldChange)
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, # plot genes with no signicant changes firstly
                                           numeric.column = 1,  # numeric.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                           cex = 0.4, pch = 16,
                                           col = add_transparency(hutils::Switch(value[[2]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                                 Not_changed="gray80",
                                                                                 Downregulated_significantly="springgreen3", 
                                                                                 Upregulated_significantly="lightcoral",
                                                                                 DEFAULT = "white"), 
                                                                  transparency = ifelse(value[[2]] == "Not_changed", 0, 1))) #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency. check https://cran.r-project.org/web/packages/circlize/circlize.pdf
                      circos.genomicPoints(region, value,
                                           numeric.column = 1,  # numeric.column=1 means the first column of value, that is the fourth of original data (het mean TPM).
                                           cex = 0.4, pch = 16,
                                           col = add_transparency(hutils::Switch(value[[2]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                                 Downregulated_significantly="springgreen3", 
                                                                                 Upregulated_significantly="lightcoral",
                                                                                 DEFAULT = "white"), 
                                                                  transparency = ifelse(value[[2]] == "Not_changed", 1, 0))) #add_transparency to add transparency of color. 0 is no transparency, 1 is total transparency. check https://cran.r-project.org/web/packages/circlize/circlize.pdf
                      circos.lines(CELL_META$cell.xlim, c(0, 0), # plot the line at y=0
                                   lty = 1, 
                                   lwd = 0.25,
                                   col = "black")
                      circos.lines(CELL_META$cell.xlim, c(-8, -8), # plot the line at min of ylim
                                   lty = 3, # lty for line type, 1 is solid line
                                   lwd = 0.3, # lwd for line thickness
                                   col = "gray80")
                      circos.lines(CELL_META$cell.xlim, c(10, 10), # plot the line at max of ylim
                                   lty = 3, 
                                   lwd = 0.3,
                                   col = "gray80")
                    })

# Step 11: Circos plot to show OR clusters to make lab log (unrelated to this project) ------------------

#####################################################################################################################
# Plot Step 11: plot Circos plot to show OR and TAAR clusters to make lab logo (not for this project)
#####################################################################################################################

#########################################################################################
################################### Start Circos plotting #############################
# not using scale_color_brewer as there are 63 different OR and TAAR clusters, but most color palettes have 8-12 colors. check http://novyden.blogspot.com/2013/09/how-to-expand-color-palette-with-ggplot.html  and https://gist.github.com/grigory93/ba4dca9636b4a6228ce5a8d5c0167968
col_for_OR_clusters_universal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(63)
names(col_for_OR_clusters_universal) <- rownames(OR_clusters_data_for_Circos_step10)
# to show the different colors with points
plot(NULL, xlim = c(1, length(col_for_OR_clusters_universal)), ylim = c(0, 4), axes = FALSE, ann = FALSE)
points(1:length(names(col_for_OR_clusters_universal)), rep(1, length(names(col_for_OR_clusters_universal))), 
       pch = 16, cex = 1, 
       col = col_for_OR_clusters_universal)
text(1:length(names(col_for_OR_clusters_universal)), rep(2.5, length(names(col_for_OR_clusters_universal))),
     names(col_for_OR_clusters_universal),
     cex=0.4, srt = 90)

# (1) plot OR clusters with chr names
pdf("/Users/qian_li/agan/Li_lab/lab_setup/Lab_recruitment/Lab_introduction/实验室介绍/Lab_logo_design/OR_TAAR_clusters_in_genome.pdf", paper='a4')
# change the basic settings of sectors and tracks
OR_and_TAAR_clusters_for_Circos_list_step11 = list(OR_clusters_data_for_Circos_step10, 
                                                   TAAR_cluster_data_for_Circos_step10)
col_for_OR_clusters_step11 <- col_for_OR_clusters_universal
circos.par(start.degree=90, 
           gap.after=0, # gap.after to adjust the gap between different sectors
           track.margin=c(0.008,0.008), # set the bottom and topmargin between tracks
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(res_all_dataframe_for_Circos_step10, plotType = "labels")
res_all_dataframe_for_Circos_step11 <- res_all_dataframe_for_Circos_step10
circos.initializeWithIdeogram(res_all_dataframe_for_Circos_step11,
                              plotType = c("labels"))  # plotType = c("axis","labels") to show axis and labels
# 1. add the first track of the different chromosomes with different colors
col_for_chr_step11 <- col_for_chr_universal
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step11,
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# 2. add the second track of position of OR and TAAR clusters as background bars
circos.genomicTrack(OR_and_TAAR_clusters_for_Circos_list_step11,
                    bg.lwd = 0.1, # adjust the thickness of border between sectors
                    bg.lty = 2, # line type, check https://r-charts.com/base-r/line-types/
                    bg.border="gray95",
                    track.height = 0.05, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      i = getI(...) # i is the index for the current data frame from the list. getI(...) has to be used inside the panel.fun.
                      sector_id = CELL_META$sector.index
                      col_for_each_region = col_for_OR_clusters_step11[match(sector_id, OR_and_TAAR_clusters_for_Circos_list_step11[[1]]$chr_name)]
                      if (i == 1) {
                        circos.genomicRect(region, value,
                                           ytop = 1, ybottom = 0,
                                           border = col_for_each_region,
                                           col = col_for_each_region)
                      } else {
                        circos.genomicRect(region, value,
                                           ytop = 1, ybottom = 0,
                                           border = "black",
                                           col = "black")
                      }
                    })
# add title
title(main = "OR and TAAR clusters (TAAR as black)",
      font.main=2, cex.main=0.8, line = -0.5) #cex.main means the magnification for main title. line means the distance between title and plot border
circos.clear()
dev.off()

# (2) only plot chr without name
pdf("/Users/qian_li/agan/Li_lab/lab_setup/Lab_recruitment/Lab_introduction/实验室介绍/Lab_logo_design/whole_genome.pdf", paper='a4')
# change the basic settings of sectors and tracks
circos.par(start.degree=90, 
           gap.after=0, # gap.after to adjust the gap between different sectors
           track.margin=c(0.008,0.008), # set the bottom and topmargin between tracks
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(res_all_dataframe_for_Circos_step10, plotType = "labels")
circos.initializeWithIdeogram(res_all_dataframe_for_Circos_step11,
                              plotType = c("labels"))  # plotType = c("axis","labels") to show axis and labels
# 1. add the first track of the different chromosomes with different colors
col_for_chr_step11 <- col_for_chr_universal
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step11,
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1
# add title
title(main = "Whole chromosome",
      font.main=2, cex.main=0.8, line = -0.5) #cex.main means the magnification for main title. line means the distance between title and plot border
circos.clear()
dev.off()


# Step 12: Circos plot to show Class I and Class II ORs and TAARs (for Yalei Kong) ------------------

#####################################################################################################################
# Plot Step 12: plot Circos plot to show Class I and Class II ORs and TAARs (not for this project)
#####################################################################################################################

#########################################################################################
################################### Start Circos plotting #############################
# (1) plot ORs and TAARs only with Class information
pdf(paste(baseDir, "/ORs_TAARs_Class_info_Circos_plot.pdf", sep=''), paper = 'a4')
# change the basic settings of sectors and tracks
circos.par(start.degree=90, 
           gap.after=0, # gap.after to adjust the gap between different sectors
           track.margin=c(0.008,0.008), # set the bottom and topmargin between tracks
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(res_all_dataframe_for_Circos_step10, plotType = "labels")
res_all_dataframe_for_Circos_step12 <- res_all_dataframe_for_Circos_step10
circos.initializeWithIdeogram(res_all_dataframe_for_Circos_step12,
                              plotType = c("labels"))  # plotType = c("axis","labels") to show axis and labels
# 1. add the first track of the different chromosomes with different colors
col_for_chr_step12 <- col_for_chr_universal
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step12,
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1

# 2. add the second track of position of Class I OR genes (blue2), class II ORs (tan) and 14 funcitonal TAAR genes (red2)
res_functional_OR_and_TAAR_step12 <- res_functional_OR_and_TAAR_step3.3 #has the column above_1Mb_between_consective_ORs
res_functional_OR_and_TAAR_dataframe_step12 <- as.data.frame(res_functional_OR_and_TAAR_step12)
res_functional_OR_and_TAAR_dataframe_for_Circos_step12 <- res_functional_OR_and_TAAR_dataframe_step12[, c("chr_name", "genomic_start", "genomic_end", "OR_Class", "zone_1to5_by_Qian_rounded_in_R")]
# add one column to indicate TAAR or OR
res_functional_OR_and_TAAR_dataframe_for_Circos_step12$is_TAAR_or_OR <- "OR"
res_functional_OR_and_TAAR_dataframe_for_Circos_step12[grepl("Taar", rownames(res_functional_OR_and_TAAR_dataframe_for_Circos_step11)), ]$is_TAAR_or_OR <- "TAAR"
# add one column to combine OR/TAAR info and class infor
res_functional_OR_and_TAAR_dataframe_for_Circos_step12$Class_info_with_OR_or_TAAR <- paste(res_functional_OR_and_TAAR_dataframe_for_Circos_step12$OR_Class, res_functional_OR_and_TAAR_dataframe_for_Circos_step12$is_TAAR_or_OR, sep = "_")
# another variable for pseudogene OR
res_pseudogene_OR_step12 <- res_pseudogene_OR_step3.2 #has the column above_1Mb_between_consective_ORs
res_pseudogene_OR_dataframe_step12 <- as.data.frame(res_pseudogene_OR_step12)
res_pseudogene_OR_dataframe_for_Circos_step12 <- res_pseudogene_OR_dataframe_step12[, c("chr_name", "genomic_start", "genomic_end", "OR_Class", "zone_1to5_by_Qian_rounded_in_R")]
# changr NA in the column of zone to 0, otherwise it may cause problems if some sectors only have NAs.
res_functional_OR_and_TAAR_dataframe_for_Circos_step12[is.na(res_functional_OR_and_TAAR_dataframe_for_Circos_step12$zone_1to5_by_Qian_rounded_in_R), ]$zone_1to5_by_Qian_rounded_in_R <- "unknown"
res_pseudogene_OR_dataframe_for_Circos_step12[is.na(res_pseudogene_OR_dataframe_for_Circos_step12$zone_1to5_by_Qian_rounded_in_R), ]$zone_1to5_by_Qian_rounded_in_R <- "unknown"
# make a list
ORs_and_TAARs_dataframe_list_for_Circos_step12 <- list(res_pseudogene_OR_dataframe_for_Circos_step12,
                                                       res_functional_OR_and_TAAR_dataframe_for_Circos_step12)
# start plotting
circos.genomicTrack(ORs_and_TAARs_dataframe_list_for_Circos_step12,
                    bg.lwd = 0.1, # adjust the thickness of border between sectors
                    bg.lty = 2,
                    bg.border="gray95",
                    track.height = 0.4, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      i = getI(...) # i is the index for the current data frame from the list. getI(...) has to be used inside the panel.fun.
                      if (i == 1) {
                        circos.genomicRect(region, value,
                                           ytop = 1, ybottom = 0,
                                           border = ifelse(value[[1]] == "Class_II", "tan", "blue2"),
                                           col = ifelse(value[[1]] == "Class_II", "tan", "blue2"))
                      } else {
                        circos.genomicRect(region, value,
                                           ytop = 1, ybottom = 0,
                                           border = hutils::Switch(value[[4]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                   Class_I_OR="blue2",
                                                                   Class_II_OR="tan", 
                                                                   NA_TAAR="red2",
                                                                   DEFAULT = "white"),
                                           col = hutils::Switch(value[[4]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                Class_I_OR="blue2",
                                                                Class_II_OR="tan", 
                                                                NA_TAAR="red2",
                                                                DEFAULT = "white"))
                      }
                    })
# 3. add the third track to show different chromosomes with different colors again
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step12,
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1

# add title
title(main = "ORs and TAARs with Class information",
      font.main=2, cex.main=0.8, line = -0.5) #cex.main means the magnification for main title. line means the distance between title and plot border
# add text to indicate gene number information
total_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)), ])
ClassI_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                          grepl("^Class_I$", res_all$OR_Class), ])
ClassI_functional_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                     grepl("^Class_I$", res_all$OR_Class) & 
                                                     grepl("protein_coding", res_all$gene_type), ])
ClassI_pseudogene_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                     grepl("^Class_I$", res_all$OR_Class) & 
                                                     grepl("pseudogene", res_all$gene_type), ])
ClassII_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                          grepl("^Class_II$", res_all$OR_Class), ])
ClassII_functional_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                     grepl("^Class_II$", res_all$OR_Class) & 
                                                     grepl("protein_coding", res_all$gene_type), ])
ClassII_pseudogene_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                     grepl("^Class_II$", res_all$OR_Class) & 
                                                     grepl("pseudogene", res_all$gene_type), ])
text(0, 0, 
     paste(total_OR_number_step12, " OR genes\n", 
           ClassI_OR_number_step12, " Class I OR genes (", 
           ClassI_functional_OR_number_step12, " functional, ", 
           ClassI_pseudogene_OR_number_step12, " pseudogenes)\n", 
           ClassII_OR_number_step12, " Class II OR genes (", 
           ClassII_functional_OR_number_step12, " functional, ", 
           ClassII_pseudogene_OR_number_step12, " pseudogenes)\n", 
           "14 functional olfactory TAAR genes, 1 TAAR psudogene", 
           sep = ""),
           cex = 0.7)
# add legend
# check Circos mannual (https://jokergoo.github.io/circlize_book/book/index.html) and ComplexHeatmap mannual (https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html)
legend_boxs_for_Class_step12 <- ComplexHeatmap::Legend(labels = c("Class I ORs", "Class II ORs", "TAARs"), 
                                                       type = "boxplot", 
                                                       legend_gp = gpar(col=c("blue2", "tan", "red2")), 
                                                       title_position = "topleft", 
                                                       title = "Track 1: Class info")
# "npc" means Normalised Parent Coordinates (the default for unit). The origin of the viewport is (0, 0) and the viewport has a width and height of 1 unit. For example, (0.5, 0.5) is the centre of the viewport.
# check https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/unit
ComplexHeatmap::draw(legend_boxs_for_Class_step12, x = unit(0.05, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
circos.clear()
dev.off()

# (2) plot ORs and TAARs with Class and zone information
pdf(paste(baseDir, "/ORs_TAARs_Class_and_zone_info_Circos_plot.pdf", sep=''), paper = 'a4')
# change the basic settings of sectors and tracks
circos.par(start.degree=90, 
           gap.after=0, # gap.after to adjust the gap between different sectors
           track.margin=c(0.008,0.008), # set the bottom and topmargin between tracks
           cell.padding = c(0,0,0,0)) #cell.padding controls the space of plot and sector border from the bottom, left, top and right 
# You can also use this one, but contain "chr" before chr names. circos.genomicInitialize(res_all_dataframe_for_Circos_step10, plotType = "labels")
circos.initializeWithIdeogram(res_all_dataframe_for_Circos_step12,
                              plotType = c("labels"))  # plotType = c("axis","labels") to show axis and labels
# 1. add the first track of the different chromosomes with different colors
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step12,
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1

# 2. add the second track of position of Class I OR genes (blue2), class II ORs (tan) and 14 funcitonal TAAR genes (red2)
circos.genomicTrack(ORs_and_TAARs_dataframe_list_for_Circos_step12,
                    bg.lwd = 0.1, # adjust the thickness of border between sectors
                    bg.lty = 2,
                    bg.border="gray95",
                    track.height = 0.2, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      i = getI(...) # i is the index for the current data frame from the list. getI(...) has to be used inside the panel.fun.
                      if (i == 1) {
                        circos.genomicRect(region, value,
                                           ytop = 1, ybottom = 0,
                                           border = ifelse(value[[1]] == "Class_II", "tan", "blue2"),
                                           col = ifelse(value[[1]] == "Class_II", "tan", "blue2"))
                      } else {
                        circos.genomicRect(region, value,
                                           ytop = 1, ybottom = 0,
                                           border = hutils::Switch(value[[4]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                   Class_I_OR="blue2",
                                                                   Class_II_OR="tan", 
                                                                   NA_TAAR="red2",
                                                                   DEFAULT = "white"),
                                           col = hutils::Switch(value[[4]], #use [[]] instead of [] to obtain vector other than dataframe.
                                                                Class_I_OR="blue2",
                                                                Class_II_OR="tan", 
                                                                NA_TAAR="red2",
                                                                DEFAULT = "white"))
                      }
                    })

# 3. add the third track of zone information
# find 5 first good colors from RColorBrewer "Accent"
col_for_OR_zones_universal <- RColorBrewer::brewer.pal(5, "Accent")
col_for_OR_zones_step12 <- col_for_OR_zones_universal
circos.genomicTrack(ORs_and_TAARs_dataframe_list_for_Circos_step12,
                    bg.lwd = 0.1, # adjust the thickness of border between sectors
                    bg.lty = 2,
                    bg.border="gray95",
                    track.height = 0.2, # set the track height, which is is the percentage to the radius of the unit circles
                    ylim = c(0,1),
                    panel.fun = function(region, value, ...) {
                      circos.genomicRect(region, value,
                                         ytop = 1, ybottom = 0,
                                         border = hutils::Switch(as.character(value[[2]]), #use [[]] instead of [] to obtain vector other than dataframe.
                                                                 "1"=col_for_OR_zones_step12[1], # hutils::Switch does not work for numeric numbers, so convert to character
                                                                 "2"=col_for_OR_zones_step12[2], 
                                                                 "3"=col_for_OR_zones_step12[3],
                                                                 "4"=col_for_OR_zones_step12[4],
                                                                 "5"=col_for_OR_zones_step12[5],
                                                                 "unknown"="NA",
                                                                 DEFAULT = "NA"),
                                         col = hutils::Switch(as.character(value[[2]]), #use [[]] instead of [] to obtain vector other than dataframe.
                                                              "1"=col_for_OR_zones_step12[1],
                                                              "2"=col_for_OR_zones_step12[2], 
                                                              "3"=col_for_OR_zones_step12[3],
                                                              "4"=col_for_OR_zones_step12[4],
                                                              "5"=col_for_OR_zones_step12[5],
                                                              "unknown"="NA",
                                                              DEFAULT = "NA"))
                    })
# Note that in the above code, you have to list "unknown"="NA", otherwise for chr12 that only have 1 pseudogene it will return nothing after hutils::Switch() function. 
# It is supposed to return DEFAULT. Maybe this is the bug for only one element that is not matched. Normally it is fine if you have multiple elements than 1.

# 4. add the fourth track to show different chromosomes with different colors again
circos.track(ylim = c(0, 1), 
             bg.col = col_for_chr_step12,
             bg.lwd = 0.01, # adjus the thickness of border
             bg.border = "gray80", 
             track.height = 0.05) #track.height is the relative height to the radius of 1

# add title
title(main = "ORs and TAARs with Class and zone information",
      font.main=2, cex.main=0.8, line = -0.5) #cex.main means the magnification for main title. line means the distance between title and plot border
# add information of gene number with different class and zone
zone1_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                         grepl("1", res_all$zone_1to5_by_Qian_rounded_in_R), ])
zone1_functional_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("1", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("protein_coding", res_all$gene_type), ])
zone1_pseudogene_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                         grepl("1", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                         grepl("pseudogene", res_all$gene_type), ])
zone2_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                         grepl("2", res_all$zone_1to5_by_Qian_rounded_in_R), ])
zone2_functional_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("2", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("protein_coding", res_all$gene_type), ])
zone2_pseudogene_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("2", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("pseudogene", res_all$gene_type), ])
zone3_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                         grepl("3", res_all$zone_1to5_by_Qian_rounded_in_R), ])
zone3_functional_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("3", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("protein_coding", res_all$gene_type), ])
zone3_pseudogene_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("3", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("pseudogene", res_all$gene_type), ])
zone4_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                         grepl("4", res_all$zone_1to5_by_Qian_rounded_in_R), ])
zone4_functional_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("4", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("protein_coding", res_all$gene_type), ])
zone4_pseudogene_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("4", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("pseudogene", res_all$gene_type), ])
zone5_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                         grepl("5", res_all$zone_1to5_by_Qian_rounded_in_R), ])
zone5_functional_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("5", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("protein_coding", res_all$gene_type), ])
zone5_pseudogene_OR_number_step12 <- nrow(res_all[grepl("^Olfr", rownames(res_all)) & 
                                                    grepl("5", res_all$zone_1to5_by_Qian_rounded_in_R) & 
                                                    grepl("pseudogene", res_all$gene_type), ])
text(0, 0, 
     paste(total_OR_number_step12, " OR genes\n", 
           ClassI_OR_number_step12, " Class I OR genes (", 
           ClassI_functional_OR_number_step12, " functional, ", 
           ClassI_pseudogene_OR_number_step12, " pseudogenes)\n", 
           ClassII_OR_number_step12, " Class II OR genes (", 
           ClassII_functional_OR_number_step12, " functional, ", 
           ClassII_pseudogene_OR_number_step12, " pseudogenes)\n", 
           "14 functional olfactory TAAR genes, 1 TAAR psudogene\n", 
           zone1_OR_number_step12, " zone 1 OR genes (", 
           zone1_functional_OR_number_step12, " functional, ", 
           zone1_pseudogene_OR_number_step12, " pseudogenes)\n", 
           zone2_OR_number_step12, " zone 2 OR genes (", 
           zone2_functional_OR_number_step12, " functional, ", 
           zone2_pseudogene_OR_number_step12, " pseudogenes)\n", 
           zone3_OR_number_step12, " zone 3 OR genes (", 
           zone3_functional_OR_number_step12, " functional, ", 
           zone3_pseudogene_OR_number_step12, " pseudogenes)\n", 
           zone4_OR_number_step12, " zone 4 OR genes (", 
           zone4_functional_OR_number_step12, " functional, ", 
           zone4_pseudogene_OR_number_step12, " pseudogenes)\n", 
           zone5_OR_number_step12, " zone 5 OR genes (", 
           zone5_functional_OR_number_step12, " functional, ", 
           zone5_pseudogene_OR_number_step12, " pseudogenes)\n", 
           sep = ""),
     cex = 0.7)
# add legend
# check Circos mannual (https://jokergoo.github.io/circlize_book/book/index.html) and ComplexHeatmap mannual (https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html)
legend_boxs_for_Class_step12 <- ComplexHeatmap::Legend(labels = c("Class I ORs", "Class II ORs", "TAARs"), 
                                                       type = "boxplot", 
                                                       legend_gp = gpar(col=c("blue2", "tan", "red2")), 
                                                       title_position = "topleft", 
                                                       title = "Track 1: Class info")
legend_boxs_for_zone_step12 <- ComplexHeatmap::Legend(labels = c("Zone 1 ORs", "Zone 2 ORs", "Zone 3 ORs", "Zone 4 ORs", "Zone 5 ORs"), 
                                                      type = "boxplot", 
                                                      legend_gp = gpar(col=col_for_OR_zones_step12), 
                                                      title_position = "topleft", 
                                                      title = "Track 2: Zone info")
# "npc" means Normalised Parent Coordinates (the default for unit). The origin of the viewport is (0, 0) and the viewport has a width and height of 1 unit. For example, (0.5, 0.5) is the centre of the viewport.
# check https://www.rdocumentation.org/packages/grid/versions/3.6.2/topics/unit
ComplexHeatmap::draw(legend_boxs_for_Class_step12, x = unit(0.05, "npc"), y = unit(0.9, "npc"), just = c("left", "center"))
ComplexHeatmap::draw(legend_boxs_for_zone_step12, x = unit(0.05, "npc"), y = unit(0.1, "npc"), just = c("left", "center"))
circos.clear()
dev.off()

# Step 13: plot counts on Trim66 exons in het and homo --------------------

#####################################################################################################################
# Plot Step 13: plot counts on Trim66 exons in het and homo
#####################################################################################################################
# import the counts on exons generated by featureCounts (see RNA_seq_code.sh for details)
counts_on_exons_featureCounts_with_gene_id_data <- read.delim(file.path(dataDir, "all_trimmomatic_featureCounts_on_exons_with_header.txt"), 
                                                              header=T, sep="\t", row.names= NULL, as.is=T)
# add gene names as a new column
counts_on_exons_featureCounts_with_gene_id_good_data <- counts_on_exons_featureCounts_with_gene_id_data[counts_on_exons_featureCounts_with_gene_id_data$Ensembl_gene_ID %in% ensembl_gene_list_good$gene_id_Ensembl, ]
counts_on_exons_featureCounts_with_gene_id_good_data$gene_name <- ensembl_gene_list_good[match(counts_on_exons_featureCounts_with_gene_id_good_data$Ensembl_gene_ID, ensembl_gene_list_good$gene_id_Ensembl), ]$gene_name    
# check if it is right
ensembl_gene_list_good[ensembl_gene_list_good$gene_name == "Trim66", ]
counts_on_exons_featureCounts_with_gene_id_good_data[counts_on_exons_featureCounts_with_gene_id_good_data$gene_name == "Trim66", ]
# extract the data on Trim66 exons
counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_step13 <- counts_on_exons_featureCounts_with_gene_id_good_data[counts_on_exons_featureCounts_with_gene_id_good_data$gene_name == "Trim66", ]
dim(counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_step13)
# import the transcript and exon information obtained by "grep Trim66 /n/data2/hms/cellbio/liberles/data/qianlong/Reference/Mus_musculus/mm10/Annotation_gtf_files/gencode.vM25.Comprehensive_gene_annotation_ALL.chr_patch_hapl_scaff.annotation.gtf > Trim66_full_info_from_Gencode_gtf_vM25.txt"
# and further processing following the linux code I used when extracting gene information from Gencode gtf file vM25 in RNA_seq_code.sh
Trim66_exon_information <- read.delim(file.path(dataDir, "Trim66_full_info_from_Gencode_gtf_vM25_with_header.txt"),
                                      header=T, sep="\t", row.names= NULL, as.is=T)
head(Trim66_exon_information)
unique(Trim66_exon_information$transcript_name) #four transcripts: "Trim66-202" "Trim66-203" "Trim66-201" "Trim66-204"
View(Trim66_exon_information) # inspect the data, it looks like that "Trim66-202" "Trim66-203" and "Trim66-201" are 3 isoforms consistent with "NM_XXX" in NCBI
Trim66_exon_information_extracted_step13 <- Trim66_exon_information[Trim66_exon_information$feature_type == "exon" & grepl("Trim66-20[1-3]", Trim66_exon_information$transcript_name), ]
# each has 20 exons, which is correct.
dim(Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-201", ])
dim(Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-202", ])
dim(Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-203", ])
# from previous research, we know that the three known isoforms (not including the two isoforms specifically expressed in MOE) have the same exons 2-18 and 20. There are slight different in exon 19.
Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-201" & Trim66_exon_information_extracted_step13$exon_number == 19, ]
Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-202" & Trim66_exon_information_extracted_step13$exon_number == 19, ]
Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-203" & Trim66_exon_information_extracted_step13$exon_number == 19, ]
# it looks like there are still differences in exon 20. The conclusion of same exon 20 is driven from comparasion of 3 NCBI isoforms, maybe it is silight different in Ensembl.
Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-201" & Trim66_exon_information_extracted_step13$exon_number == 20, ]
Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-202" & Trim66_exon_information_extracted_step13$exon_number == 20, ]
Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-203" & Trim66_exon_information_extracted_step13$exon_number == 20, ]
# So I will use Trim66-201 (NM_181853.4) as an model to extract exon 2-20 for downstream analysis
Trim66_201_exon_information_extracted_step13 <- Trim66_exon_information_extracted_step13[Trim66_exon_information_extracted_step13$transcript_name == "Trim66-201", ]
dim(Trim66_201_exon_information_extracted_step13)
# does the exons match with the featureCounts data? The answer is yes.
all(Trim66_201_exon_information_extracted_step13$genomic_start %in% counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_step13$genomic_start & 
      Trim66_201_exon_information_extracted_step13$genomic_end %in% counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_step13$genomic_end)
which(counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_step13$genomic_start %in% Trim66_201_exon_information_extracted_step13$genomic_start &
        counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_step13$genomic_end %in% Trim66_201_exon_information_extracted_step13$genomic_end)
# use dplyr to merge the data together for exons of Trim66-201. There are same rows (in terms of exons) with same exon fragment counts. This is becuase they belong to different transcripts and featureCounts still list them separately, but with same amount of fragment counts.
counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_201_isoform_step13 <- counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_step13 %>% 
  dplyr::distinct() %>% # keep the unique rows
  left_join(Trim66_201_exon_information_extracted_step13) %>% # obtain the full exon information
  filter(!is.na(exon_id)) %>% #remove rows with NAs
  arrange(exon_number)

# melt the data for downstream analysis using gather
counts_on_exons_Trim66_201_isoform_melted_grouped_step13 <- counts_on_exons_featureCounts_with_gene_id_good_data_Trim66_201_isoform_step13 %>% 
  dplyr::select(chr_name:genomic_end, Trim66_KO_het_1:Trim66_KO_homo_5, exon_number:exon_id) %>% 
  tidyr::gather(key ="samples", value="fragment_counts", Trim66_KO_het_1:Trim66_KO_homo_5) %>% # gather to melt data for plot
  mutate(genotype = ifelse(grepl("het", samples), "het", "homo"))

# obtain the summary data
counts_on_exons_Trim66_201_isoform_melted_grouped_summary_step13 <- counts_on_exons_Trim66_201_isoform_melted_grouped_step13 %>% 
  group_by(exon_number, genotype) %>% 
  summarize(mean_counts=mean(fragment_counts), 
            sd_counts=sd(fragment_counts), 
            Numbers=n(), 
            se=sd_counts/sqrt(mean_counts), 
            upper_limit=mean_counts+se, 
            lower_limit=mean_counts-se)

# delete exon 1 data as exon 1 is not the same between different isoforms
counts_on_exons_Trim66_201_isoform_melted_grouped_summary_delete_exon1_step13 <- 
  counts_on_exons_Trim66_201_isoform_melted_grouped_summary_step13 %>% 
  filter(exon_number != 1)
counts_on_exons_Trim66_201_isoform_melted_grouped_delete_exon1_step13 <- 
  counts_on_exons_Trim66_201_isoform_melted_grouped_step13 %>% 
  filter(exon_number != 1)

# obtain the unpaired t test results
ttest_results_Trim66_exons_2to20 <- data.frame(matrix(nrow = 19, ncol = 2))
colnames(ttest_results_Trim66_exons_2to20) <- c("exon_number", "unpaired_t_test")
ttest_results_Trim66_exons_2to20$exon_number <- paste("exon No.", 2:20, sep = "")
for (i in 2:20) {
  ttest_results_temp <- t.test(
    counts_on_exons_Trim66_201_isoform_melted_grouped_step13[grepl("het", counts_on_exons_Trim66_201_isoform_melted_grouped_step13$samples) & 
                                                               counts_on_exons_Trim66_201_isoform_melted_grouped_step13$exon_number == i, ]$fragment_counts, 
    counts_on_exons_Trim66_201_isoform_melted_grouped_step13[grepl("homo", counts_on_exons_Trim66_201_isoform_melted_grouped_step13$samples) & 
                                                               counts_on_exons_Trim66_201_isoform_melted_grouped_step13$exon_number == i, ]$fragment_counts,
    alternative = "two.sided",
    var.equal = TRUE) #use var.equal = TRUE, otherwise, the values are larger
  ttest_results_Trim66_exons_2to20[i-1, "unpaired_t_test"] <- ttest_results_temp$p.value
}

# start plotting
plot_step13 <- ggplot() + 
  geom_bar(counts_on_exons_Trim66_201_isoform_melted_grouped_summary_delete_exon1_step13,
           mapping=aes(x=exon_number, y=mean_counts, color=genotype),
           stat="identity", 
           position = position_dodge(0.95),
           width = 0.9, #adjust the bar width
           size=0.25, # adjust the line thickness
           fill="gray95") + 
  geom_errorbar(counts_on_exons_Trim66_201_isoform_melted_grouped_summary_delete_exon1_step13,
                mapping=aes(x=exon_number, ymin=lower_limit, ymax=upper_limit, color=genotype),
                position = position_dodge(0.95), width = .2, size=0.25) +
  scale_color_manual(values = c("springgreen3", "red1")) +
  geom_point(counts_on_exons_Trim66_201_isoform_melted_grouped_delete_exon1_step13,
             mapping=aes(x=exon_number, y=fragment_counts, group=genotype, fill=genotype),
             position=position_jitterdodge(jitter.width = 0.2,
                                           dodge.width = 0.95),
             size=3, shape=21, stroke = 0.4, # stroke control border thickness
             color="black", alpha =0.5) +
  scale_fill_manual(values = rep("white", 2)) +
  theme(legend.position = "top",
        plot.title = element_text(size = rel(2)), # adjust the title size
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.5))) +
  ggtitle('Fragment counts on Trim66 exons (2-20, shared exons among 3 known NM_XXX isoforms) in Trim66 KO (5 samples each) by ggplot') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("Exon number (from Trim66-201 transcript Ensembl ID: ENSMUST00000033339, NCBI ID: NM_181853)") + 
  ylab("Raw fragment counts (mean +- se)") +
  scale_x_continuous(breaks = 2:20) + # label 2 to 20 other than the default 5, 10, 20
  annotate("text", x = 2, y = 700, 
           label = paste(capture.output(as.data.frame(ttest_results_Trim66_exons_2to20)), sep=" ",collapse = "\n"),
           size=4, colour = "black") + # if the printed text have several lines, adjust the console size, as capture.output function is similar to take the screenshot of in console
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

plot_step13
ggsave(paste(baseDir, "/results/exon_counts_in_Trim66_KO.pdf", sep=''), plot_step13, width = 50, height = 25, units ="cm")

#####################################################################################################################
                # Summary of data processing tips: 
########  Always use tidyverse package and pipeline!!!!
########  tidyverse package makes life much easier !!!
#####################################################################################################################


# Step 14: plot receptor changes against Ngn1+/Omp+ or Gap43+/Gap43- -------------------------------

#####################################################################################################################
# Plot Step 14: plot OR and TAAR changes against RPKM FC of Ngn1+/Omp+ (Stavtos data 2011, Cell, Lsd1)  
# or TPM FC of Gap43+/Gap43- (our data, fixed cells stained with Gap43 then FACS) to see if they are correlated
# The idea is to check if ORs expressed early in neuron progenitors of immature OSNs will be upregulated
#####################################################################################################################
library(readxl)
# import the original TPM and RPKM data, note that some genes with name March.. may be messed. Not an issue here.
TPM_Gap43_cells_data_Longzhi <- readxl::read_excel("/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/151215_TPM_data_sorted_by_Longzhi.xlsx", 
                                                   sheet = 1)
RPKM_Ngn1_Omp_cells_data_Stavros <- readxl::read_excel("/Users/qian_li/agan/Li_lab/Results/Qian_data/Deep_sequencing/RNA_seq/Database/Olfactory_receptors/151215_TPM_data_sorted_by_Longzhi.xlsx", 
                                                       sheet = 3)
# change the column names, change the gene names (change -ps to .ps to be consistent with res_all rownames), 
# and add a column as FC of Gap43+/Gap43- or Ngn1+/Omp+ (add 0.1 to Gap43- and Omp+ to avoid 0s)
TPM_Gap43_cells_data_Longzhi_step14 <- TPM_Gap43_cells_data_Longzhi %>% 
  rename(Gap43_positive = `Gap43+`, Gap43_negative = `Gap43-`) %>% 
  dplyr::select(gene:Gap43_negative) %>% 
  distinct(gene, .keep_all=TRUE) %>% # remove dupliacted rows with the same gene name
  mutate(gene_name_converted = str_replace(gene, "-", "."), FC_Gap43_positive_to_negative = Gap43_positive/(Gap43_negative + 0.1)) 
RPKM_Ngn1_Omp_cells_data_Stavros_step14 <- RPKM_Ngn1_Omp_cells_data_Stavros %>% 
  rename(Ngn1_positive = `Ngn1+`, Omp_positive = `Omp+`) %>% 
  dplyr::select(gene:Omp_positive) %>% 
  distinct(gene, .keep_all=TRUE) %>% 
  mutate(gene_name_converted = str_replace(gene, "-", "."), FC_Ngn1_to_Omp = Ngn1_positive/(Omp_positive + 0.1))
# How many functional ORs are in the above 2 dataset? And add columns with log2FC and significance
TPM_Gap43_cells_data_Longzhi_functional_ORs_step14 <- TPM_Gap43_cells_data_Longzhi_step14 %>% 
  filter(gene_name_converted %in% rownames(res_functional_OR_and_TAAR_dataframe_step6.1)) %>% 
  left_join(dplyr::select(res_functional_OR_and_TAAR_dataframe_step6.1, log2FoldChange, padj, if_significant_upregulate_or_downregulate, name), 
            by=c("gene_name_converted" = "name")) %>% #add the FC and padj from res
  mutate(log2FoldChange_NA_to_0 = ifelse(is.na(log2FoldChange), 0, log2FoldChange), 
         padj_NA_to_1 = ifelse(is.na(padj), 1, padj)) %>% #change NA in log2FC to 0 and NA in padj to 1
  filter(grepl("Olfr", gene_name_converted)) # only select functional ORs
dim(TPM_Gap43_cells_data_Longzhi_functional_ORs_step14)
sum(rownames(res_functional_OR) %in% TPM_Gap43_cells_data_Longzhi_functional_ORs_step14$gene_name_converted)

RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14 <- RPKM_Ngn1_Omp_cells_data_Stavros_step14 %>% 
  filter(gene_name_converted %in% rownames(res_functional_OR_and_TAAR_dataframe_step6.1)) %>% 
  left_join(dplyr::select(as.data.frame(res_functional_OR_and_TAAR_dataframe_step6.1), log2FoldChange, padj, if_significant_upregulate_or_downregulate, name), 
            by=c("gene_name_converted" = "name")) %>% #add the FC and padj from res
  mutate(log2FoldChange_NA_to_0 = ifelse(is.na(log2FoldChange), 0, log2FoldChange), 
         padj_NA_to_1 = ifelse(is.na(padj), 1, padj)) %>% #change NA in log2FC to 0 and NA in padj to 1
  filter(grepl("Olfr", gene_name_converted)) # only select functional ORs
dim(RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14)
sum(rownames(res_functional_OR) %in% RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14$gene_name_converted)

# how many TAARs are in the above 2 dataset
TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14 <- TPM_Gap43_cells_data_Longzhi_step14 %>% 
  filter(gene_name_converted %in% rownames(res_functional_OR_and_TAAR_dataframe_step6.1)) %>% 
  left_join(dplyr::select(res_functional_OR_and_TAAR_dataframe_step6.1, log2FoldChange, padj, if_significant_upregulate_or_downregulate, name), 
            by=c("gene_name_converted" = "name")) %>% #add the FC and padj from res
  mutate(log2FoldChange_NA_to_0 = ifelse(is.na(log2FoldChange), 0, log2FoldChange), 
         padj_NA_to_1 = ifelse(is.na(padj), 1, padj)) %>% #change NA in log2FC to 0 and NA in padj to 1
  filter(grepl("Taar", gene_name_converted)) # only select functional TAARs
dim(TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14)

RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14 <- RPKM_Ngn1_Omp_cells_data_Stavros_step14 %>% 
  filter(gene_name_converted %in% rownames(res_functional_OR_and_TAAR_dataframe_step6.1)) %>% 
  left_join(dplyr::select(as.data.frame(res_functional_OR_and_TAAR_dataframe_step6.1), log2FoldChange, padj, if_significant_upregulate_or_downregulate, name), 
            by=c("gene_name_converted" = "name")) %>% #add the FC and padj from res
  mutate(log2FoldChange_NA_to_0 = ifelse(is.na(log2FoldChange), 0, log2FoldChange), 
         padj_NA_to_1 = ifelse(is.na(padj), 1, padj)) %>% #change NA in log2FC to 0 and NA in padj to 1
  filter(grepl("Taar", gene_name_converted)) # only select functional TAARs
dim(RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14)

# How many pseudogene ORs are in the above 2 dataset
TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14 <- TPM_Gap43_cells_data_Longzhi_step14 %>% 
  filter(gene_name_converted %in% rownames(res_pseudogene_OR_dataframe_step6.2)) %>% 
  left_join(dplyr::select(res_pseudogene_OR_dataframe_step6.2, log2FoldChange, padj, if_significant_upregulate_or_downregulate, name), 
            by=c("gene_name_converted" = "name")) %>% #add the FC and padj from res
  mutate(log2FoldChange_NA_to_0 = ifelse(is.na(log2FoldChange), 0, log2FoldChange), 
         padj_NA_to_1 = ifelse(is.na(padj), 1, padj)) #change NA in log2FC to 0 and NA in padj to 1
dim(TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14)
sum(rownames(res_pseudogene_OR) %in% TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14$gene_name_converted)

RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14 <- RPKM_Ngn1_Omp_cells_data_Stavros_step14 %>% 
  filter(gene_name_converted %in% rownames(res_pseudogene_OR_dataframe_step6.2)) %>% 
  left_join(dplyr::select(as.data.frame(res_pseudogene_OR_dataframe_step6.2), log2FoldChange, padj, if_significant_upregulate_or_downregulate, name), 
            by=c("gene_name_converted" = "name")) %>% #add the FC and padj from res
  mutate(log2FoldChange_NA_to_0 = ifelse(is.na(log2FoldChange), 0, log2FoldChange), 
         padj_NA_to_1 = ifelse(is.na(padj), 1, padj)) #change NA in log2FC to 0 and NA in padj to 1
dim(RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14)
sum(rownames(res_pseudogene_OR) %in% RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14$gene_name_converted)

# 14.1: plot the receptor changes against Gap43+/Gap43-
pdf(paste(resultDir, "/ORs_TAARs_changes_Trim66_KO_against_gap43+_to_gap43-_FC.pdf", sep=''), width = 100/2.54, height = 100/2.54)
# functional ORs
plot_step14.1_1 <- ggplot(data=TPM_Gap43_cells_data_Longzhi_functional_ORs_step14, aes(x=FC_Gap43_positive_to_negative, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 50, y = 5, 
           label = paste("Pearson Correlation: ", 
                         round(cor(TPM_Gap43_cells_data_Longzhi_functional_ORs_step14$log2FoldChange_NA_to_0, TPM_Gap43_cells_data_Longzhi_functional_ORs_step14$FC_Gap43_positive_to_negative), 10),
                         ", p value: ",
                         round(cor.test(TPM_Gap43_cells_data_Longzhi_functional_ORs_step14$log2FoldChange_NA_to_0, TPM_Gap43_cells_data_Longzhi_functional_ORs_step14$FC_Gap43_positive_to_negative)$p.value, 10), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 agains Gap43+/Gap43-') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("FC of receptors in Gap43+ vs Gap43-") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
plot_step14.1_2 <-  plot_step14.1_1+
  geom_text(data=subset(TPM_Gap43_cells_data_Longzhi_functional_ORs_step14, !(if_significant_upregulate_or_downregulate == "Not_changed")), 
            aes(label = gene_name_converted), 
            size =5, colour='blue')  # label=name in the aes() suggest the text to label

arrange_step14.1_1and2 <- ggarrange(plot_step14.1_1, plot_step14.1_2, ncol = 1, nrow = 2)

# functional TAARs
plot_step14.1_3 <- ggplot(data=TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14, aes(x=FC_Gap43_positive_to_negative, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 3, y = 1, 
           label = paste("Pearson Correlation: ", 
                         round(cor(TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14$log2FoldChange_NA_to_0, TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14$FC_Gap43_positive_to_negative), 10),
                         ", p value: ",
                         round(cor.test(TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14$log2FoldChange_NA_to_0, TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14$FC_Gap43_positive_to_negative)$p.value, 10), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 agains Gap43+/Gap43-') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("FC of receptors in Gap43+ vs Gap43-") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
plot_step14.1_4 <-  plot_step14.1_3+
  geom_text(data=subset(TPM_Gap43_cells_data_Longzhi_functional_TAARs_step14, !(if_significant_upregulate_or_downregulate == "Not_changed")), 
            aes(label = gene_name_converted), 
            size =5, colour='blue')  # label=name in the aes() suggest the text to label

arrange_step14.1_3and4 <- ggarrange(plot_step14.1_3, plot_step14.1_4, ncol = 1, nrow = 2)

# pseudogene ORs
plot_step14.1_5 <- ggplot(data=TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14, aes(x=FC_Gap43_positive_to_negative, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 20, y = 1, 
           label = paste("Pearson Correlation: ", 
                         round(cor(TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14$log2FoldChange_NA_to_0, TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14$FC_Gap43_positive_to_negative), 10),
                         ", p value: ",
                         round(cor.test(TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14$log2FoldChange_NA_to_0, TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14$FC_Gap43_positive_to_negative)$p.value, 10), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 agains Gap43+/Gap43-') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("FC of receptors in Gap43+ vs Gap43-") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
plot_step14.1_6 <-  plot_step14.1_5+
  geom_text(data=subset(TPM_Gap43_cells_data_Longzhi_pseudogene_ORs_step14, !(if_significant_upregulate_or_downregulate == "Not_changed")), 
            aes(label = gene_name_converted), 
            size =5, colour='blue')  # label=name in the aes() suggest the text to label

arrange_step14.1_5and6 <- ggarrange(plot_step14.1_5, plot_step14.1_6, ncol = 1, nrow = 2)

arrange_step14.1_1and2
arrange_step14.1_3and4
arrange_step14.1_5and6
dev.off()

# 14.2: plot the receptor changes against Ngn1+/Omp+
pdf(paste(resultDir, "/ORs_TAARs_changes_Trim66_KO_against_Ngn1+_to_Omp+_FC.pdf", sep=''), width = 100/2.54, height = 100/2.54)
# functional ORs
plot_step14.2_1 <- ggplot(data=RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14, aes(x=FC_Ngn1_to_Omp, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 6, y = 5, 
           label = paste("Pearson Correlation: ", 
                         round(cor(RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14$log2FoldChange_NA_to_0, RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14$FC_Ngn1_to_Omp), 10),
                         ", p value: ",
                         round(cor.test(RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14$log2FoldChange_NA_to_0, RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14$FC_Ngn1_to_Omp)$p.value, 10), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 agains Ngn1+/Omp+') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("FC of receptors in Ngn1+ vs Omp+") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
plot_step14.2_2 <-  plot_step14.2_1+
  geom_text(data=subset(RPKM_Ngn1_Omp_cells_data_Stavros_functional_ORs_step14, !(if_significant_upregulate_or_downregulate == "Not_changed")), 
            aes(label = gene_name_converted), 
            size =5, colour='blue')  # label=name in the aes() suggest the text to label

arrange_step14.2_1and2 <- ggarrange(plot_step14.2_1, plot_step14.2_2, ncol = 1, nrow = 2)

# functional TAARs
plot_step14.2_3 <- ggplot(data=RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14, aes(x=FC_Ngn1_to_Omp, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 0.1, y = 1, 
           label = paste("Pearson Correlation: ", 
                         round(cor(RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14$log2FoldChange_NA_to_0, RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14$FC_Ngn1_to_Omp), 10),
                         ", p value: ",
                         round(cor.test(RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14$log2FoldChange_NA_to_0, RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14$FC_Ngn1_to_Omp)$p.value, 10), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Functional TAARs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 agains Ngn1+/Omp+') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("FC of receptors in Ngn1+ vs Omp+") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
plot_step14.2_4 <-  plot_step14.2_3+
  geom_text(data=subset(RPKM_Ngn1_Omp_cells_data_Stavros_functional_TAARs_step14, !(if_significant_upregulate_or_downregulate == "Not_changed")), 
            aes(label = gene_name_converted), 
            size =5, colour='blue')  # label=name in the aes() suggest the text to label

arrange_step14.2_3and4 <- ggarrange(plot_step14.2_3, plot_step14.2_4, ncol = 1, nrow = 2)

# pseudogene ORs
plot_step14.2_5 <- ggplot(data=RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14, aes(x=FC_Ngn1_to_Omp, y=log2FoldChange_NA_to_0)) + 
  geom_point(aes(fill=if_significant_upregulate_or_downregulate, #interaction to group two columns in ggplot2. check https://stackoverflow.com/questions/9968976/group-by-two-columns-in-ggplot2
                 color=if_significant_upregulate_or_downregulate,
                 alpha=if_significant_upregulate_or_downregulate), # if the legend does not have filled color, also check https://github.com/tidyverse/ggplot2/issues/2322
             size=8, shape=21) + #The graphical argument used to specify point shapes is "pch", 21 is filled circle so we can change the filled color. http://www.sthda.com/english/wiki/r-plot-pch-symbols-the-different-point-shapes-available-in-r
  stat_smooth(formula = y ~ x, method = lm, colour = "blue") + # add Pearson correlation ±95% confidence interval
  scale_color_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_fill_manual(values = c("springgreen3", "gray40", "red1")) +
  scale_alpha_manual(values = c(1,0.5,1)) +
  geom_hline(yintercept=0, linetype="solid", 
             color = "black", size=1) +
  annotate("text", x = 4, y = 1, 
           label = paste("Pearson Correlation: ", 
                         round(cor(RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14$log2FoldChange_NA_to_0, RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14$FC_Ngn1_to_Omp), 10),
                         ", p value: ",
                         round(cor.test(RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14$log2FoldChange_NA_to_0, RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14$FC_Ngn1_to_Omp)$p.value, 10), sep = ""),
           hjust = 1.1, vjust = -.5, size=10, colour = "blue", fontface =2) + #add the text about the pearson correlation coeficient
  theme(legend.position = "top",
        plot.title = element_text(size = rel(3)), # adjust the title size
        axis.title = element_text(size = rel(2)),
        axis.text = element_text(size = rel(2))) +
  ggtitle('Pseudogene ORs changed in Trim66 KO homo vs het with threhold padj<0.05 and abs(FoldChange)>=1.5 agains Ngn1+/Omp+') +
  theme(plot.title = element_text(vjust = 1)) + # adjust the title position
  xlab("FC of receptors in Ngn1+ vs Omp+") + ylab("log2 gene expression changes of homo to het") +
  labs(color = "homo vs het significant", fill = "homo vs het significant", alpha = "homo vs het significant") + # change legend title
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) # remove background grey color and border lines et al. check http://felixfan.github.io/ggplot2-remove-grid-background-margin/

# add text
plot_step14.2_6 <-  plot_step14.2_5+
  geom_text(data=subset(RPKM_Ngn1_Omp_cells_data_Stavros_pseudogene_ORs_step14, !(if_significant_upregulate_or_downregulate == "Not_changed")), 
            aes(label = gene_name_converted), 
            size =5, colour='blue')  # label=name in the aes() suggest the text to label

arrange_step14.2_5and6 <- ggarrange(plot_step14.2_5, plot_step14.2_6, ncol = 1, nrow = 2)

arrange_step14.2_1and2
arrange_step14.2_3and4
arrange_step14.2_5and6
dev.off()

