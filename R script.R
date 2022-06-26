library(Seurat)
library(ggplot2)
library(AnnotationHub)
library(SingleCellExperiment)
library(pheatmap)
library(limma)
library(edgeR)


setwd("H:/Apply/Data Regeneration/Practice SCE/6/")
dex <- readRDS("rds/dex.rds")
V.Hours_Dex <- c("00","01","02","04","08","18")


# Seurat Object Creation ####
# Reading the raw data that downloaded from the GEO database
dex <- read.delim("Data/GSE141834_scRNAseq_rawCounts.txt", header = T, row.names = 1, sep = "\t")
# Creating a seurat object from the UMI count matrix
dex <- CreateSeuratObject(dex,project = "DEX")
# Normalizing the data by log10 normalization method and scaling by 10000
dex <- NormalizeData(dex, normalization.method = "LogNormalize", scale.factor = 10000)
# Adding the time course of Dexamethasone treatment for each cell
dex$Hours_Dex <- substr(colnames(dex), 5, 6)
# Calculating the percentage of mitochondrial genes in each cell
dex$percent.mt <- PercentageFeatureSet(dex, pattern = "^MT-")
# Removing the cells with the percentage of mitochondrial genes higher than 5%
# as said in the paper
dex <- subset(dex, subset = percent.mt < 5)


# Plotting ####
# Fig. 1a
ggplot(mapping = aes(dex$Hours_Dex,dex$nCount_RNA)) + 
  geom_violin(aes(fill = factor(dex$Hours_Dex))) + 
  geom_jitter(size= 0.3)+
  theme_classic()+
  theme(legend.title = element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15, face = "bold"))+
  labs(x="Hours Dex", y="Unique Transcripts/Cell")
# Fig. 1b
ggplot(mapping = aes(dex$Hours_Dex,dex$nFeature_RNA)) + 
  geom_violin(aes(fill = factor(dex$Hours_Dex))) + 
  geom_jitter(size= 0.3)+
  theme_classic()+
  theme(legend.title = element_blank(), 
        axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 15, face = "bold"))+
  labs(x="Hours Dex", y="Features/Cell")



# Calculating the ratio of cells that expressing special genes 
# in each time course of Dex treatment
genes <- c("PNMT","TSC22D3","GLUL","ZNF703","SPTSSB", "IER3", "PLK2")
exp.ratio <- data.frame()
x <- NULL
for (j in genes) {
  for (i in 1:length(V.Hours_Dex)) {
    x <- dex[j]
    x <- x[,c(dex$Hours_Dex == V.Hours_Dex[i] & x$nCount_RNA != 0)]
    exp.ratio <-rbind(exp.ratio , c(j , V.Hours_Dex[i],dim(x)[2]/400))
  }
}
colnames(exp.ratio) <- c("gene", "Hours_Dex", "Ratio")
#saveRDS(exp.ratio, "rds/exp.ratio.rds")


# Scaled Log10 expression (box and scatter plots) and
# ratio of expressing cells (black dot/line plot)
# Fig. 1f
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')["PNMT",], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')["PNMT",], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_line(inherit.aes = F, mapping = aes(x = exp.ratio[exp.ratio$gene == "PNMT",]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == "PNMT",]$Ratio)*3, group =1)) +
  geom_point(mapping = aes(x = exp.ratio[exp.ratio$gene == "PNMT",]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == "PNMT",]$Ratio)*3)) +
  scale_y_continuous(sec.axis = sec_axis(~ ./3, breaks = seq(0,1,0.25),
                                         name = "Ratio of Expressing Cells")) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = "PNMT") +
  theme_classic() + NoLegend()
# Fig. 1g
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[2],], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[2],], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_line(inherit.aes = F, mapping = aes(x = exp.ratio[exp.ratio$gene == genes[2],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[2],]$Ratio)*3, group =1)) +
  geom_point(mapping = aes(x = exp.ratio[exp.ratio$gene == genes[2],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[2],]$Ratio)*3)) +
  scale_y_continuous(sec.axis = sec_axis(~ ./3, breaks = seq(0,1.25,0.25),
                                         name = "Ratio of Expressing Cells")) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = genes[2]) +
  theme_classic() + NoLegend()
# Fig. 1h
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[3],], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[3],], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_line(inherit.aes = F, mapping = aes(x = exp.ratio[exp.ratio$gene == genes[3],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[3],]$Ratio)*4, group =1)) +
  geom_point(mapping = aes(x = exp.ratio[exp.ratio$gene == genes[3],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[3],]$Ratio)*4)) +
  scale_y_continuous(sec.axis = sec_axis(~ ./4, breaks = seq(0,1,0.25),
                                         name = "Ratio of Expressing Cells")) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = genes[3]) +
  theme_classic() + NoLegend()
# Fig. 2e
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[4],], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[4],], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_line(inherit.aes = F, mapping = aes(x = exp.ratio[exp.ratio$gene == genes[4],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[4],]$Ratio)*3, group =1)) +
  geom_point(mapping = aes(x = exp.ratio[exp.ratio$gene == genes[4],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[4],]$Ratio)*3)) +
  scale_y_continuous(sec.axis = sec_axis(~ ./3, breaks = seq(0,1,0.25),
                                         name = "Ratio of Expressing Cells")) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = genes[4]) +
  theme_classic() + NoLegend()
# Supplemental Fig. 1b
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[5],], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[5],], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_line(inherit.aes = F, mapping = aes(x = exp.ratio[exp.ratio$gene == genes[5],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[5],]$Ratio)*2, group =1)) +
  geom_point(mapping = aes(x = exp.ratio[exp.ratio$gene == genes[5],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[5],]$Ratio)*2)) +
  scale_y_continuous(sec.axis = sec_axis(~ ./2, breaks = seq(0,1,0.25),
                                         name = "Ratio of Expressing Cells")) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = genes[5]) +
  theme_classic() + NoLegend()
# Supplemental Fig. 1c
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[6],], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[6],], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_line(inherit.aes = F, mapping = aes(x = exp.ratio[exp.ratio$gene == genes[6],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[6],]$Ratio)*2.75, group =1)) +
  geom_point(mapping = aes(x = exp.ratio[exp.ratio$gene == genes[6],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[6],]$Ratio)*2.75)) +
  scale_y_continuous(sec.axis = sec_axis(~ ./2.75, breaks = seq(0,1,0.25),
                                         name = "Ratio of Expressing Cells")) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = genes[6]) +
  theme_classic() + NoLegend()
# Supplemental Fig. 1d
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[7],], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')[genes[7],], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_line(inherit.aes = F, mapping = aes(x = exp.ratio[exp.ratio$gene == genes[7],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[7],]$Ratio)*3, group =1)) +
  geom_point(mapping = aes(x = exp.ratio[exp.ratio$gene == genes[7],]$Hours_Dex, y = as.numeric(exp.ratio[exp.ratio$gene == genes[7],]$Ratio)*3)) +
  scale_y_continuous(sec.axis = sec_axis(~ ./3, breaks = seq(0,1,0.25),
                                         name = "Ratio of Expressing Cells")) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = genes[7]) +
  theme_classic() + NoLegend()


# Cell Cycle Scoring ####

# AnnotationHub
ah <- AnnotationHub()
# Access the Ensembl database for organism
ahDb <- query(ah, pattern = c("Homo sapiens", "EnsDb"), ignore.case = TRUE)
# Acquire the latest annotation files
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
# Download the appropriate Ensembldb database
edb <- ah[[id]]
# Extract gene-level information from database
annotations <- genes(edb, return.type = "data.frame")
# Select annotations of interest
ens_id <- annotations %>% dplyr::select(gene_id, gene_name)
#saveRDS(ens_id, "ens_id.rds")


# Converting HS.Pairs' ENSEMBL ID to gene symbol
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))

rn.dex <- data.frame(g1f = hs.pairs$G1$first, g1s = hs.pairs$G1$second)
rn.dex <- dplyr::left_join(rn.dex, ens_id, by = c("g1f" = "gene_id"))
rn.dex <- dplyr::left_join(rn.dex, ens_id, by = c("g1s" = "gene_id"))
colnames(rn.dex)[3:4] <- c("first", "second")
hs.pairs$G1 <- rn.dex[3:4]

rn.dex <- data.frame(sf = hs.pairs$S$first, ss = hs.pairs$S$second)
rn.dex <- dplyr::left_join(rn.dex, ens_id, by = c("sf" = "gene_id"))
rn.dex <- dplyr::left_join(rn.dex, ens_id, by = c("ss" = "gene_id"))
colnames(rn.dex)[3:4] <- c("first", "second")
hs.pairs$S <- rn.dex[3:4]

rn.dex <- data.frame(g2mf = hs.pairs$G2M$first, g2ms = hs.pairs$G2M$second)
rn.dex <- dplyr::left_join(rn.dex, ens_id, by = c("g2mf" = "gene_id"))
rn.dex <- dplyr::left_join(rn.dex, ens_id, by = c("g2ms" = "gene_id"))
colnames(rn.dex)[3:4] <- c("first", "second")
hs.pairs$G2M <- rn.dex[3:4]


# Substituting the old gene symbols of the Seurat object with the new ones
new.rownames <- data.frame(rn = rownames(dex))
new.rownames <- dplyr::left_join(new.rownames, ens_id, by = c("rn" = "gene_name"))
new.rownames <- new.rownames[!duplicated(new.rownames$rn),]

na.genes <- new.rownames[is.na(new.rownames$gene_id) | grepl("LRG",new.rownames$gene_id),]
na.genes$subs <- limma::alias2SymbolTable(na.genes$rn, species = 'Hs')
na.genes <- na.genes[!is.na(na.genes$subs),]
na.genes <- dplyr::left_join(na.genes, ens_id, by = c("subs" = "gene_name"))
na.genes <- na.genes[!grepl("LRG", na.genes$gene_id.y),]
na.genes <- na.genes[!duplicated(na.genes$subs),]

new.rownames <- dplyr::left_join(new.rownames, na.genes, by = c("rn" = "rn"))
new.rownames[!is.na(new.rownames$subs),]$rn <- new.rownames[!is.na(new.rownames$subs),]$subs
new.rownames[!is.na(new.rownames$gene_id.y),]$gene_id <- new.rownames[!is.na(new.rownames$gene_id.y),]$gene_id.y
new.rownames <- new.rownames[,-3:-5]
new.rownames <- dplyr::left_join(new.rownames, ens_id, by = c("gene_id" = "gene_id"))
new.rownames[is.na(new.rownames$gene_name),]$gene_name <- new.rownames[is.na(new.rownames$gene_name),]$rn
new.rownames <- new.rownames[,-1]
#saveRDS(new.rownames, "new.rownames.rds")


# Converting Seurat object into a Singlecellexperiment object
dex.sce <- as.SingleCellExperiment(dex)
# Cell cycle scoring by means of the Scran package based on the converted row names
cycle.scores <- scran::cyclone(dex.sce, hs.pairs, gene.names = new.rownames$gene_name)


# Transferring the assigned cell cycle scores into the Seurat object
dex$G1 <- cycle.scores$scores$G1
dex$S <- cycle.scores$scores$S
dex$G2M <- cycle.scores$scores$G2M


# Scaling the Seurat object ####
# Finding the top 500 variable genes for scaling
dex <- FindVariableFeatures(dex, nfeatures = 500)
# Scaling and centering the top 500 variable genes using the negative binomial model
# The impacts of transcript counts, mitochondrial percent, and cell-cycle scores 
# has been regressed-out.
dex <- ScaleData(dex, model.use = "negbinom",
                 vars.to.regress = 
                   c("percent.mt", 
                     "G1","S","G2M", 
                     "nCount_RNA"))


# Principal Components, tSNE and UMAP ####
set.seed(1234)
# Performing principal component analysis using the top 30 principal components
dex <- RunPCA(dex, verbose = F, npcs = 30, seed.use = 119)
DimPlot(dex, reduction = "pca", group.by = "Hours_Dex", pt.size = 1.5) +
  scale_y_continuous(breaks = seq(-15,10,5))+
  theme(legend.title = element_text(size = 20),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17)) + 
  guides(color= guide_legend(title = "Hrs Dex"))

# Performing tSNE using the top 16 principal components
dex <- RunTSNE(dex, reduction = "pca", dims = 1:16, perplexity = 40)
# Fig. 3c
DimPlot(dex, reduction = "tsne", group.by = "Hours_Dex", pt.size = 1.5) +
  scale_y_continuous(breaks = seq(-20,30,10))+
  scale_x_continuous(breaks = seq(-30,20,10))+
  theme(legend.title = element_text(size = 20),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17)) + 
  guides(color= guide_legend(title = "Hrs Dex"))
# Note: A slight difference between this plot and the t-SNE plot in the paper is due to
# the fact that we don't know what value has been set to perplexity and seed.use in 
# the RunTSNE function

# Performing UMAP using the top 16 principal components
dex <- RunUMAP(dex, reduction = "pca", dims = 1:16)
# Supplemental Fig. 4a
DimPlot(dex, reduction = "umap", group.by = "Hours_Dex", pt.size = 1.5) +
  theme(legend.title = element_text(size = 20),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17))
# Note: A slight difference between this plot and the UMAP plot in the paper is due to
# the fact that we don't know what value has been set to the parameters used in 
# the RunUMAP function



# Feature plots ####
featureplot.genes <- c("FKBP5","PNMT","TSC22D3","HSD11B2","ZNF703","DDIT4")
# Supplemental Fig. 6a
FeaturePlot(object = dex, features = featureplot.genes, pt.size = 0.9, 
            reduction = "tsne", slot = "data", ncol = 3)


# Clustering ####
# Finding nearest neighbors using the top 16 principal components
dex <- FindNeighbors(dex, reduction = "pca", dims = 1:16)
# Clustering the data into 7 clusters using a proper resolution
dex <- FindClusters(dex, resolution = 0.45)
#saveRDS(dex, "rds/dex.rds")


# Calculating the number of cells from each treatment timepoint present in the clusters
cell.clust <- data.frame()
cell.clust <- hist_fun()
# Fig. 3g
ggplot(cell.clust, aes(Hours_Dex, as.numeric(value), fill =Hours_Dex)) + 
  geom_bar(stat='identity') + facet_wrap(~variable, ncol = 7)+
  scale_y_continuous(breaks = c(0,50,100,150,200,250,300))


hist_fun <- function(){
  a <- data.frame()
  c <- NULL
  b <- NULL
  a <- rownames(1:6)
  for (i in 0:6) {
    b <- dex[,dex$seurat_clusters == i]$Hours_Dex
    for (j in 1:length(V.Hours_Dex)) {
      c[j] <- length(b[b==V.Hours_Dex[j]])
    }
    a <- cbind(a , as.numeric(c))
  }
  print(1)
  a <- cbind(a, V.Hours_Dex)
  colnames(a) <- c(paste0("cluster",0:6), "Hours_Dex")
  a <- as.data.frame(a)
  a <- a[,c(6,3,1,2,7,5,4,8)]
  a <- reshape2::melt(a, id.vars = "Hours_Dex")
  return(a)
}



# Differentially Expressed Genes ####
# Finding the differentially expressed genes by comparing each dex-treatment timepoint
# with the untreated cells. 
# For this matter, we exploit the MAST test with a fold change cutoff of 1.25, 
# an adjusted p-value of 0.01, and excluding genes detected in fewer than 10% of cells.
table.DEGs <- data.frame()
sc.DEGs <- NULL
sc.DEGs.table <- NULL
for (i in 2:6) {
  sc.DEGs.table[[paste0("dex",V.Hours_Dex[i])]] <- 
    FindMarkers(dex, ident.1 = V.Hours_Dex[i], ident.2 = "00", 
                group.by = "Hours_Dex",test.use = "MAST", slot = "data")
  b <- sc.DEGs.table[[paste0("dex",V.Hours_Dex[i])]]
  b$padj <- p.adjust(b$p_val, method = "BH")
  sc.DEGs.table[[paste0("dex",V.Hours_Dex[i])]] <- 
    b[b$padj < 0.01 & b$pct.1 > 0.1 & abs(b$avg_logFC) > log10(1.25),]
  sc.DEGs <- append(sc.DEGs, rownames(sc.DEGs.table[[paste0("dex",V.Hours_Dex[i])]]))
  up <- NULL
  up <- dim(b[b$padj < 0.01 & b$pct.1 > 0.1 & b$avg_logFC > log10(1.25),])[1]
  down <- NULL
  down <- dim(b[b$padj < 0.01 & b$pct.1 > 0.1 & b$avg_logFC < -log10(1.25),])[1]
  table.DEGs <- rbind(table.DEGs, c(up+down,up,down))
}
colnames(table.DEGs) <- c("# DEGs", "# up","# down")
rownames(table.DEGs) <- c("01","02","04","08","18")
table.DEGs$singlecell <- t(table.DEGs$singlecell)
# Receiving the unique list of DEGs across all timepoints
DEGs <- NULL
DEGs$singlecell <- unique(sc.DEGs)
# Calling single-cell DEGs which not exist in bulk DEGs
DEGs$singlecell.uniq <- DEGs$singlecell[!(DEGs$singlecell %in% DEGs$bulk)]



# Bulk RNA Sequencing ####
# reading the normalized data as a count matrix
cn <- read.delim("Data/GSE141834_bulkRNAseq_normalized_counts.txt", header = T,
                 row.names = 1, sep = "\t")
# Un-normalize the data
cn <- (2^(cn))


## Bulk-RNA Differential Expression Analysis ####
# Constructing the design matrix for Differential Expression Analysis
design <- NULL
design <- cbind("dex00"=c(rep(1,3),rep(0,15)),
                "dex01"=c(rep(0,3),rep(1,3),rep(0,12)),
                "dex02"=c(rep(0,6),rep(1,3),rep(0,9)),
                "dex04"=c(rep(0,9),rep(1,3),rep(0,6)),
                "dex08"=c(rep(0,12),rep(1,3),rep(0,3)),
                "dex18"=c(rep(0,15),rep(1,3)))
rownames(design) <- colnames(cn)
# Performing differential expression analysis by means of Limma-voom package
v <- voom(counts = round(cn)-1, design =design, normalize = "quantile")
fit <- lmFit(v, design)
cont.matrix <- NULL
cont.matrix$bulk01 <- makeContrasts(dex01-dex00, levels = design)
cont.matrix$bulk02 <- makeContrasts(dex02-dex00, levels = design)
cont.matrix$bulk04 <- makeContrasts(dex04-dex00, levels = design)
cont.matrix$bulk08 <- makeContrasts(dex08-dex00, levels = design)
cont.matrix$bulk18 <- makeContrasts(dex18-dex00, levels = design)

# Calling and filtering DEGs for each Dex-treatment timepoint with a fold change 
# cutoff of 1.5 and adjusted p-value of 0.05.
table.bulk.DEGs <- data.frame()
bulk.DEGs.df <- list()
bulk.DEGs <- NULL
fig2c.heatmap.df <- NULL
for (i in 2:6) {
  fit2 <- contrasts.fit(fit, cont.matrix[[paste0("bulk", V.Hours_Dex[i])]])
  fit2 <- eBayes(fit2)
  bulk.DEGs.df[[paste0("bulk", V.Hours_Dex[i])]] <- 
    topTable(fit2, sort.by = "B", number = Inf)
  
  b <- bulk.DEGs.df[[paste0("bulk", V.Hours_Dex[i])]]
  
  bulk.DEGs.df[[paste0("bulk", V.Hours_Dex[i])]] <- 
    b[b$adj.P.Val < 0.05 & abs(b$logFC) > log2(1.5),]
  
  #if(i == 3){
  #  a <- b[b$adj.P.Val < 0.05,]
  #  bulk.DEGs.df[[paste0("bulk", V.Hours_Dex[i])]] <- 
  #    a[a$logFC > log2(1.5) | a$logFC < -1.5,]
  #}
  
  
  fig2c.heatmap.df[[paste0("bulk", V.Hours_Dex[i])]] <-
    b[DEGs$singlecell.uniq,]
  
  
  bulk.DEGs <- append(bulk.DEGs, 
                      rownames(bulk.DEGs.df[[paste0("bulk",V.Hours_Dex[i])]]))
  up <- NULL
  up <- dim(b[b$adj.P.Val < 0.05 & b$logFC > log2(1.5),])[1]
  down <- NULL
  down <- dim(b[b$adj.P.Val < 0.05 & b$logFC < -log2(1.5),])[1]
  table.bulk.DEGs <- rbind(table.bulk.DEGs, c(up+down,up,down))
  
  
}
colnames(table.bulk.DEGs) <- c("# DEGs", "# up","# down")
rownames(table.bulk.DEGs) <- c("01","02","04","08","18")
# Adding the resulted number of DEGs to the table of DEGs
table.DEGs$bulk <- t(table.bulk.DEGs)
#saveRDS(table.DEGs, "rds/table.DEGs.rds")

# Receiving the unique list of DEGs across all timepoints
DEGs$bulk <- unique(bulk.DEGs)
# Calling bulk DEGs which not exist in single-cell DEGs
DEGs$bulk.uniq <- DEGs$bulk[!(DEGs$bulk %in% DEGs$singlecell)]
#saveRDS(DEGs, "rds/DEGs.rds")



## Plotting ####
# Plotting Fig. 1e and supplemental Fig. 1a
# Calculating mean of log10 values for each gene in each timepoint
cn.log10 <- log10(round(cn)-1)
bulk.row.mean <- data.frame("00" = rowMeans(cn.log10[1:3]),
                "01" = rowMeans(cn.log10[4:6]),
                "02" = rowMeans(cn.log10[7:9]),
                "04" = rowMeans(cn.log10[10:12]),
                "08" = rowMeans(cn.log10[13:15]),
                "18" = rowMeans(cn.log10[16:18]))

# Calculating standerd deviation for each gene at each timepoint for plotting error bar
cn.colnames <- c("EtOH", "D1hr","D2hr","D4hr","D8hr","D18hr")
bulk.plot.genes <- c("PNMT", "TSC22D3", "GLUL","SPTSSB","IER3","PLK2")
bulk.sds <- NULL
for (i in 1:length(cn.colnames)) {
  a <- cn.log10[,colnames(cn.log10)[grepl(cn.colnames[i], colnames(cn.log10))]]
  sda <- NULL
  for (j in 1:length(bulk.plot.genes)) {
    sda <- append(sda, sd(a[bulk.plot.genes[j],]))
  }
  bulk.sds <- rbind(bulk.sds, sda)
}
rownames(bulk.sds) <- V.Hours_Dex
colnames(bulk.sds) <- bulk.plot.genes


table.fig1e <- as.data.frame(t(bulk.row.mean[c("PNMT", "TSC22D3", "GLUL"),]))
table.fig1e$Hours_Dex <- c(0, 1, 2, 4, 8, 18)
table.fig1e <- reshape2::melt(table.fig1e, id.vars = "Hours_Dex")
table.fig1e$sd <- reshape2::melt(bulk.sds)$value[1:18]

table.sfig1a <- as.data.frame(t(bulk.row.mean[c("SPTSSB","IER3","PLK2"),]))
table.sfig1a$Hours_Dex <- c(0, 1, 2, 4, 8, 18)
table.sfig1a <- reshape2::melt(table.sfig1a, id.vars = "Hours_Dex")
table.sfig1a$sd <- reshape2::melt(bulk.sds)$value[19:36]

# Fig 1e
ggplot(table.fig1e, aes(Hours_Dex, value, color = variable)) +
  geom_line(size= 1) + geom_point() + theme_classic() +
  geom_errorbar(aes(ymin= value-sd ,ymax= value+sd), width=0.4, color="black", size= 0.9) +
  labs(x="Hours Dex", y="Bulk RNAseq log10(Mean Reads)") +
  guides(color = guide_legend(title = "Gene"))

# Supplemental Fig. 1a
ggplot(table.sfig1a, aes(Hours_Dex, value, color = variable)) +
  geom_line(size= 1) + geom_point() + theme_classic() +
  geom_errorbar(aes(ymin= value-sd ,ymax= value+sd), width=0.3, color="black", size= 0.9) +
  labs(x="Hours Dex", y="Bulk RNAseq log10(Mean Reads)") +
  guides(color = guide_legend(title = "Gene"))



# Plotting DEGs ####
# Scaling and centering each DEG
dex.deg <- ScaleData(dex, features = DEGs$singlecell, 
                     vars.to.regress = c("percent.mt", 
                                         "G1","S","G2M","nCount_RNA"))
# Change idents to dex-treatment timepoints
Idents(dex.deg) <- dex.deg$Hours_Dex
avg.exp <- AverageExpression(dex.deg, slot = "scale.data",
                       features = DEGs$singlecell.uniq)$RNA
hm.sc <- pheatmap(avg.exp, border_color = NA, cluster_cols = F, show_rownames = F,
                  fontsize_col = 20, angle_col = 0)
row.order <- hm.sc$tree_row$order

# Plotting Fig. 2c
fig2c.heatmap.df <- NULL
for (i in 2:6) {
  fit2 <- contrasts.fit(fit, cont.matrix[[paste0("bulk", V.Hours_Dex[i])]])
  fit2 <- eBayes(fit2)
  fig2c.heatmap.df[[paste0("bulk", V.Hours_Dex[i])]] <- 
    topTable(fit2, sort.by = "B", number = Inf)
  
  fig2c.heatmap.df[[paste0("bulk", V.Hours_Dex[i])]] <- 
    fig2c.heatmap.df[[paste0("bulk", V.Hours_Dex[i])]][DEGs$singlecell.uniq,]
  
  
}

bulk.df.heatmap <- NULL
bulk.df.heatmap <- data.frame(bulk00 = rep(0,174))
for (i in 2:6) {
  bulk.df.heatmap <- cbind(bulk.df.heatmap,  
                           fig2c.heatmap.df[[paste0("bulk", V.Hours_Dex[i])]]$logFC)
}
colnames(bulk.df.heatmap) <- V.Hours_Dex
rownames(bulk.df.heatmap) <- DEGs$singlecell.uniq
bulk.df.heatmap[is.na(bulk.df.heatmap)] <- 0
# Fig. 2c
phmb <- pheatmap(bulk.df.heatmap[row.order,], cluster_rows = F, 
                 cluster_cols = F, show_rownames = F, 
                 fontsize_col = 20, angle_col = 0)



# Ratio of Responding Genes ####

# Calculating the mean log-scaled expression level and the standard deviation
# for plotting the horizontal line of Fig. 4a and 4b
dex.sce <- as.SingleCellExperiment(dex)
dex.fkbp5 <- dex.sce["FKBP5"][,dex.sce$Hours_Dex=="00"]
# Excluding zero values form calculation
dex.fkbp5 <- dex.fkbp5[,which(counts(dex.fkbp5) != 0)]
hline.fkbp5 <- mean(logcounts(dex.fkbp5)) + sd(logcounts(dex.fkbp5))

dex.ier3 <- dex.sce["IER3"][,dex.sce$Hours_Dex=="00"]
dex.ier3 <- dex.ier3[,which(counts(dex.ier3) != 0)]
hline.ier3 <- mean(logcounts(dex.ier3)) - sd(logcounts(dex.ier3))

# Fig. 4a
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')["FKBP5",], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')["FKBP5",], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_hline(yintercept = hline.fkbp5, color = "darkblue", size = 1) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = "FKBP5") +
  theme_classic() + NoLegend()

# Fig. 4b
ggplot() +
  geom_boxplot(mapping = aes(y= GetAssayData(dex, slot = 'data')["IER3",], x= dex$Hours_Dex, 
                             colour = dex$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= GetAssayData(dex, slot = 'data')["IER3",], x= dex$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  geom_hline(yintercept = hline.ier3, color = "darkblue", size = 1) +
  labs(x= "Hours Dex", y= "Scaled Log10 Expression", title = "IER3") +
  theme_classic() + NoLegend()



# Determining how many Dex target genes showed a response in each cell (Ratio of Responding Genes)
# We determine the mean log-scaled expression level and standard deviation 
# for each Differentially Expressed Gene in untreated cells (using only non-zero values).
# A gene is "responsive" if it was expressed greater than one SD above the mean
# untreated level or more than one SD below the mean for downregulated genes.
rrg.df <- data.frame()
for (i in 1:414) {
  print(i)
  a <- NULL
  a <- dex.sce[sc.DEGs[i]][,dex.sce$Hours_Dex=="00"]
  a <- a[,which(counts(a) != 0)]
  a18 <- dex.sce[sc.DEGs[i]][,dex.sce$Hours_Dex=="18"]
  a18 <- a18[,which(counts(a18) != 0)]
  sda <- sd(logcounts(a))
  if (is.na(sda)) {
    sda = 0
  }
  mean.a <- mean(logcounts(a))
  mean.a18 <- mean(logcounts(a18))
  if (is.na(mean.a)) {
    mean.a = 0
  }
  if (is.na(mean.a18)) {
    mean.a18 = 0
  }
  cond = NULL
  if (mean.a18 < mean.a) {
    cond = "dr"
  } else {
    cond = "ur"
  }
  b <- logcounts(dex.sce[sc.DEGs[i]])
  e.list <- list()
  for (j in 1:2400) {
    if (mean.a == 0 & mean.a18 == 0) {
      e.list <- append(e.list, 3)
    } else if (b[j] == 0) {
      e.list <- append(e.list, 3)
    } else {
      if (cond == "ur") {
        if (b[j] >= mean.a + sda) {
          e.list <- append(e.list, 1)
        } else {
          e.list <- append(e.list, 0)
        }
      } else {
        if (b[j] <= mean.a - sda) {
          e.list <- append(e.list, 1)
        } else {
          e.list <- append(e.list, 0)
        }
      }
    }
  }
  rrg.df <- rbind(rrg.df, e.list)
}
#saveRDS(rrg.df, "rds/rrg.df.rds")


# Generating a list that contains the RRG for each cell
rrg.list <- NULL
for (i in 1:2400) {
  rrg.list <- append(rrg.list, 
                     length(rrg.df[i][rrg.df[i] == 1])/length(rrg.df[i][rrg.df[i] != 3]))
}
# Adding the RRG list to the Seurat object as a meta data
dex[["RRG"]] <- rrg.list

# Fig 4c
ggplot() +
  geom_boxplot(mapping = aes(y= dex$RRG, x= dex.sce$Hours_Dex, 
                             colour = dex.sce$Hours_Dex),outlier.shape = NA) +
  geom_jitter(mapping = aes(y= dex$RRG, x= dex.sce$Hours_Dex),
              size = 0.2, alpha = 0.3, width = 0.30) +
  labs(x= "Hours Dex", y= "Ratio of Responding Genes") +
  theme_classic() + NoLegend()

# Fig 4d
FeaturePlot(dex, features = "RRG", reduction = "pca", pt.size = 1.3)+
  scale_y_continuous(breaks = seq(-15,10,5))+
  scale_color_gradientn(colors = c("#0011ff", "#0088ff","#00FFFF", "#00FF40", 
                                   "#80FF00", "#FFFF00", "#FFBF00", "#FF0000", 
                                   "#FF00BF", "#FF00BF")) + 
  theme(legend.title = element_text(size = 20), 
        plot.title = element_blank(),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17))

# Fig 4e
FeaturePlot(dex, features = "RRG", reduction = "tsne", pt.size = 1.3)+
  scale_y_continuous(breaks = seq(-20,30,10))+
  scale_x_continuous(breaks = seq(-30,20,10))+
  scale_color_gradientn(colors = c("#0011ff", "#0088ff","#00FFFF", "#00FF40", 
                                   "#80FF00", "#FFFF00", "#FFBF00", "#FF0000", 
                                   "#FF00BF", "#FF00BF")) + 
  theme(legend.title = element_text(size = 20), 
        plot.title = element_blank(),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17))
