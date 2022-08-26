# Load required packages

#Cran
library("ggplot2")
library("RColorBrewer")
library("circlize")
# Bioconductor
library("DESeq2")
library("apeglm")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("ReactomePA")
library("enrichplot")
#Github
library("ComplexHeatmap")
library("EnhancedVolcano")


# Set working directory (make sure to change this to where you cloned the repo)
setwd("GeneExpressionWorkshopRMed2022")

# Create a directory to store output files
dir.create("output")
dir.create("output/figures")

# Get input data
read_counts <- read.csv("data/ProstateCancer_DMSO_SP2509_LSD1i_readcounts.csv")
sampleInfo <- read.csv("data/ProstateCancer_sampleInfo.csv")

# Inspect the data
View(read_counts)
View(sampleInfo)

# Set read counts rownames to gene names and drop gene_id column (we don't need it)
rownames(read_counts) <- read_counts$gene_id
read_counts <- read_counts[,-1]

# To eliminate errors with zero read counts, we increment all reads by 1
read_counts <- read_counts + 1
head(read_counts)

# Set the sample info rownames to sample names
rownames(sampleInfo) <- sampleInfo$sample
head(sampleInfo)

# Some checks to make sure read counts table and sample info table are properly prepared
check1 <- all(rownames(sampleInfo) %in% colnames(read_counts))
check2 <- all(rownames(sampleInfo) == colnames(read_counts))
temp_counts <- read_counts[, rownames(sampleInfo)]
check3 <- all(rownames(sampleInfo) == colnames(temp_counts))

# These should all return true if everything was prepared correctly
check1
check2
check3



### Create DESeq2 data set (dds) ###

# We specify the names of the conditions/groups to used as the reference and the treatment
control_factor <- "DMSO"
treatment_factor <- "SP2509"

# Construct a DESeqDataSet (dds) using the read count table, sample info table, and specify the reference condition
dds <- DESeqDataSetFromMatrix(countData = read_counts, 
                              colData = sampleInfo, 
                              design = ~ condition)

# Let's view the dds object
View(dds)

# For now, we only have one assay, which are the raw reads
assays(dds)

# The raw read counts for each sample and gene can be obtained by the counts() method
View(counts(dds))
View(assays(dds)$counts)

# Every column that was present in the sample info table (```sampleInfo```) is also built into the dds object as metadata columns.
# We can access it using the column name as shown. This will be handy later on to compare gene expression data with respect to clinical attributes.
dds$condition
dds$replicate
dds$psa

# Usually we only want to keep genes that have at least 10 reads (or some other minimum threshold amount).
# For large experiments with many samples, this can speed up analysis, but is not required.
# Let's filter the dds object to only keep genes that have at least 10 reads
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]

# Re-order the condition levels so that the reference factor ("DMSO" in our case) is first
dds$condition <- relevel(dds$condition, ref = control_factor)

# Perform differential expression analysis
dds <- DESeq(dds)

# View the dds object
dds

# Checking the DESeq Analysis by visualizing
# Dispersion Estimate Plot
plotDispEsts(dds,ylim=c(1e-6, 1e2))

# Extract the results from the dds object. This is done using the results() method
res <- results(dds, contrast=c("condition", treatment_factor, control_factor))

# View it better as a dataframe
DF_res <- data.frame(res@listData)
rownames(DF_res) <- rownames(dds)
View(DF_res)

# MA plot Visualization
plotMA(res, ylim=c(-3,3))

summary (res)

# Perform log fold change(LFC) shrinkage using ```apeglm```. This looks at large LFC that are not due to 
# low read counts and uses that to inform a prior distribution. Genes that have large LFC but high statistical info
# are not shrunk, whereas LFC with little statistical data are shrunk to reduce bias. This useful for plotting the data, 
# and also for comparing gene expression data from one independent experiment with other independent experiments.
LFC_coef <- paste0("condition_", treatment_factor, "_vs_", control_factor)
resLFC <- lfcShrink(dds, coef=LFC_coef, type="apeglm")

#Visualize Shrunk data
plotMA(resLFC, ylim=c(-3,3))


# Filter results based on a False Discovery Rate (FDR) cut off.
# Set the FDR (false discovery rate) percentage value expressed as a decimal.
# Typically this is 10% (0.1) or 5% (0.05) if you want to be more stringent.
FDR_aplha <- 0.05
resFDR <- results(dds, alpha=FDR_aplha)
head(resFDR)

# Since all of the gene IDs are currently specified as ENSEMBL IDs (typically what's provided after read alignment),
# we need to map them to HGNC gene symbols so we can easily identify the genes. 
# We use the Bioconductor libraries ```AnnotationDbi``` and ```org.Hs.eg.db``` (since this is human data) to do this.
# We will map to HGNC gene symbols and store this in a ```GeneID``` column in both ```res```  and ```resFDR```.
res$GeneID <- mapIds(org.Hs.eg.db, keys=rownames(res), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
resFDR$GeneID <- mapIds(org.Hs.eg.db, keys=rownames(resFDR), column="SYMBOL", keytype="ENSEMBL", multiVals="first")

# Let's re-order the columns so ```GeneID``` is first
res <- res[,c(7,1:6)]
resFDR <- resFDR[,c(7,1:6)]

# View the results with mapped gene symbols
head(res)
head(resFDR)

# View it better as a dataframe
DF_res <- data.frame(res@listData)
rownames(DF_res) <- rownames(dds)

# Let's make a separate table of ENSEMBL symbols mapped to gene ids. This will be useful later in downstream applications
geneList <- as.data.frame(res)
geneList <- subset(geneList, select = c(GeneID)) # Keep only the GeneID column; drop everything else
geneList$ENSEMBL <- rownames(geneList) # Save ENSEMBL IDs in new ENSEMBL column

View(geneList)

# Save results tables
write.csv(as.data.frame(res), file='output/results.csv')
write.csv(as.data.frame(resFDR), file='output/results_FDR.csv')



# For visualizing the data, it is useful to calculate a variance stabilizing transformation (VST) from the 
# fitted dispersion model. yielding a matrix of values which are now approximately homoskedastic (having constant 
# Variance along the range of mean values). We need to set ```blind``` to ```FALSE``` here to prevent re-estimation of dispersion
# and maintain any large differences in LFC/counts due to inherent biology of the experiment.
vst <- varianceStabilizingTransformation(dds, blind = FALSE)

# View the vst object
vst
View(assay(vst))




### Draw PCA plot ### 

# Use the plotPCA() method to perform PCA and get the output data
pcaData <- plotPCA(vst, intgroup=c("condition", "replicate"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Generate the plot using ggplot2
pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
            geom_point(size = 3) +
            labs(shape="Replicate", color="Condition") +
            xlab(paste0("PC1: ",percentVar[1],"% variance")) +
            ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
            theme_classic()
pca_plot

# Save to PDF
pdf("output/figures/pca_plot.pdf", width=5, height=5)
pca_plot
dev.off()




### Plot a sample clustering using ComplexHeatmap ###

# Transpose the data
sampleDists <- dist(t(assay(vst)))

# Prepare a matrix with the data
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst$condition, vst$replicate, sep="-")
colnames(sampleDistMatrix) <- paste(vst$condition, vst$replicate, sep="-")
map_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Generate heatmap for sample clustering
clustering_plot <- Heatmap(sampleDistMatrix,
                           name = "Variance",
                           col = map_colors,
                           #row distance method
                           clustering_distance_rows = "pearson",
                           #column distance method
                           clustering_distance_columns = "pearson",
                           show_row_names = TRUE,
                           show_column_names = TRUE,
                           row_names_side = "left",
                           column_names_side = "bottom")
clustering_plot

# Save to PDF
pdf("output/figures/sample_clustering_plot.pdf", width=10, height=8)
clustering_plot
dev.off()




### Draw gene expression heatmap using ComplexHeatmap ###

# Re-order the vst data from most differentially expressed to least differentially expressed.
# Let's pull out only the top 30 differentially expressed genes.
genestokeep <- order(rowVars(assay(vst)), decreasing = TRUE)[1:30]
heatmap_data <- as.data.frame(assay(vst)[genestokeep,])
head(heatmap_data)

# Get the HGNC gene symbols for the subsetted data
heatmap_data$GeneID <- mapIds(org.Hs.eg.db,keys=rownames(heatmap_data),column="SYMBOL",keytype="ENSEMBL",multiVals="first")

# Inspect the results. You will notice that some ENSEMBL IDs don't map to a HGNC gene symbol (they don't have gene names).
head(heatmap_data)

# Let's remove them by removing NAs
heatmap_data <- na.omit(heatmap_data, cols = c("GeneID"))

# Set rownames to ```GeneID```
rownames(heatmap_data) <- heatmap_data$GeneID

# We need to drop the ```GeneID``` column before sending to ComplexHeatmap because the input must be numerical
heatmap_data <- subset(heatmap_data, select = -c(GeneID))

# Scale the data across all genes so differences among samples are not overpowered by strong outliers.
heatmap_data2 = t(apply(heatmap_data, 1, function(x) {scale(x)}))
colnames(heatmap_data2) = colnames(heatmap_data)
#View(heatmap_data)
#View(heatmap_data2)
head(heatmap_data2)


# Create a HeatmapAnnotation object to specify the colours for conditions and replicates.
replicate_colors = c("firebrick1", "dodgerblue", "darkseagreen", "darkviolet", "lightsalmon")
condition_colors = c("grey", "aquamarine3")
names(replicate_colors) = unique(sampleInfo[,"replicate"])
names(condition_colors) = unique(sampleInfo[,"condition"]) 
heatmap_anno = HeatmapAnnotation(Replicates = as.matrix(colData(dds)[,c("replicate")]),
                                 Condition = as.matrix(colData(dds)[,c("condition")]), 
                                 col = list(Replicates = replicate_colors), Condition = condition_colors)


#generate heatmap
expression_heatmap <- Heatmap(heatmap_data2,
                              name = "Expression",
                              top_annotation = heatmap_anno,
                              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                              cluster_rows = TRUE,
                              cluster_columns = FALSE,
                              #row distance method
                              clustering_distance_rows = "spearman",
                              #column distance method
                              clustering_distance_columns = "spearman",
                              #row cluster method
                              clustering_method_rows = "centroid",
                              #column cluster method
                              clustering_method_columns = "centroid",
                              row_names_side = "left",
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              column_names_side = "bottom",
                              #position of row dendrogam
                              row_dend_side = "left",
                              #position of column dendrogram
                              column_dend_side = "top",
                              #width of row dendrogram
                              row_dend_width = unit(1, "cm"),
                              #height of column dendrogram
                              column_dend_height = unit(1, "cm")
                              )
expression_heatmap

# Save to PDF
pdf("output/figures/expression_heatmap.pdf", width=7, height=8)
expression_heatmap
dev.off()




### Draw Volcano plot using EnhancedVolcano ###

# Generate a volcano plot of all genes
# Here we can specify cutoffs for LFC (1.0) and p values (0.05) and colour code genes accordingly
# Genes coloured red have greater than +/-  1.0 LFC and are statistically significant.
volcano_plot <- EnhancedVolcano(resFDR,
                                title = paste(control_factor, treatment_factor, sep = " vs. "),
                                subtitle = "",
                                lab = resFDR$GeneID,
                                x = "log2FoldChange",
                                y = "padj",
                                pCutoff = 0.05,
                                FCcutoff = 1,
                                gridlines.major = FALSE,
                                gridlines.minor = FALSE,
                                legendPosition = "bottom",
                                legendLabels = c('Not Significant', expression(Log[2]~FC~only), "p-value only", expression(p-value~and~log[2]~FC)),
                                cutoffLineType = "longdash",
                                cutoffLineCol = 'black',
                                labCol = 'black',
                                col = c("#5B5B5C", "#0FBD32", "#126FE8", "#FF0000"),
                                caption = ""
                                )
volcano_plot

# Save to PDF
pdf("output/figures/volcano_plot.pdf", width=8, height=8)
volcano_plot
dev.off()




### Individual gene plots ###

# Gene to plot
gene_to_plot <- "CDKN1A"
gene_to_plot <- "VEGFA"
gene_to_plot <- "CCNB1"
gene_to_plot <- "BRCA1"
gene_to_plot <- "BRCA2"
gene_to_plot <- "EZH2"

# Because we are using the HGNC gene symbol and the dds object references genes by ENSEMBL ID, we lookup the HGNC gene name
# in our previously created ```geneList``` table to get the corresponding ENSEMBL ID. Of course, we could have also used an 
# ENSEMBL ID directly.
# Get normalized read counts from dds object for selected gene
gene_counts <- plotCounts(dds, gene=geneList[which(geneList$GeneID==gene_to_plot),2], intgroup=c("condition", "replicate", "mutationCount", "psa", "tmb"), returnData=TRUE)
head(gene_counts)

# We can plot normalized read counts and visualize individual replicates
p1 <- ggplot(gene_counts, aes(x = condition, y = count)) +
    ggtitle(gene_to_plot) +
    geom_jitter(aes(color = replicate), size = 2) +
    xlab("") +
    ylab("Normalized read count") +
    theme_classic()
p1

# We can plot normalized read counts and see how it corresponds to tumor mutational burden
p2 <- ggplot(gene_counts, aes(x = condition, y = count)) +
    ggtitle(gene_to_plot) +
    geom_jitter(aes(color = tmb), size = 2) +
    xlab("") +
    ylab("Normalized read count") +
    theme_classic()
p2

# We can plot normalized read counts and see how it corresponds to PSA levels at the time of biopsy
p3 <- ggplot(gene_counts, aes(x = condition, y = count)) +
    ggtitle(gene_to_plot) +
    geom_jitter(aes(color = psa), size = 2) +
    xlab("") +
    ylab("Normalized read count") +
    theme_classic()
p3

# Save to PDF
pdf("output/figures/normalized_expression_plots.pdf", width=5, height=5)
p1
p2
p3
dev.off()





### Perform Pathway Enrichment Analysis ###

# Before performing any enrichment analysis, we need to filter our results even further to select genes of interest.
# We can do this by subsetting genes that have a specific LFC and p value
# Let's subset genes that have a LFC greater than 2.0 and an adjusted p value smaller than 0.05.
filtered_resFDR <- as.data.frame(resFDR)

# Let's remove any genes that have ```NA``` for adjusted p value
filtered_resFDR <- na.omit(filtered_resFDR, cols = c("padj"))

# Subset the gene table
filtered_resFDR <- subset(filtered_resFDR,
                          log2FoldChange >= 2 &
                            padj < 0.05
)
head(filtered_resFDR)

# For ```enrichplot``` package requires gene names to be in ENTREZID format, so we need to convert from ENSEMBL to ENTREZID
filtered_resFDR$EntrezID <- mapIds(org.Hs.eg.db,keys=rownames(filtered_resFDR),column="ENTREZID",keytype="ENSEMBL",multiVals="first")

# Calculate pathway enrichment
enrichment_data <- enrichPathway(gene = filtered_resFDR$EntrezID, pvalueCutoff = 0.05, readable = TRUE)

enrichment_data2 <- pairwise_termsim(enrichment_data)
# Generate enrichment map
enrichment_plot <- emapplot(enrichment_data2) +
  scale_color_continuous(low = "#E62412", high = "#374AB3", trans = "reverse", labels = to_scientific) +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.text.align = 0)

# Modify legend
enrichment_plot <- enrichment_plot + guides(colour = "colorbar", size = "legend")
enrichment_plot <- enrichment_plot + guides(size = guide_legend(title = "Gene count"))
enrichment_plot <- enrichment_plot + guides(color = guide_colorbar(title = "P adjusted"))

tree_plot <- treeplot(enrichment_data2)
heat_plot <- heatplot(enrichment_data2, showCategory=5)

# Draw plot
tree_plot
heat_plot
enrichment_plot


