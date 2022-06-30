# Analysis of Single-Cell and Bulk RNA-seq Data by R
Rewriting the code of an article about employing scRNA-seq and bulk RNA-seq for research on breast cancer. 
<br> This project has been done mainly by means of the Seurat package (version 3) in R.</br>
<br> The main challenge of this project was to rewrite the code without knowing the value of parameters that have been set in various functions. As a result, I had to run a function several times to find the best value for a parameter.
## Table of Contents
- [Introduction](#introduction)
- [The Main Goal](#the-main-goal)
- [Brief Summary of The Article](#brief-summary-of-the-article)
- [Softwares and Packages](#softwares-and-packages)
- [Data Analysis and Plotting Figures](#data-analysis-and-plotting-figures)

## Introduction
Cancer is one of the leading causes of death in the world and it can affect people of all ages. Cancer cells are a result of genetic mutations caused by endogenous or environmental factors. The development of malignant tumors is a complex process that is mostly uncovered and it needs to make manifold efforts to unravel the main pathway that every normal cell passes until it turns into a cancer cell. Obviously, shedding light on this multi-stage process to decipher the hidden aspects of cancer cells' progression.
<br> Poor prognosis, the spread of cancer cells into distant organs due to metastasis, evasion from the immune system, and resistance to medications are among the major hurdles in the way of cancer therapy. For this reason, researchers in this field have to get familiar with cancer cells and tumor microenvironment at different stages in order to uncover the vulnerable points of the disease for designing brand-new medications. 
<br> Nowadays, manifold sequencing techniques and methods can make an important contribution to achieving this goal. Especially, single-cell RNA sequencing (scRNA-seq) and bulk RNA sequencing are featured methods that can help to investigate cellular heterogeneity of tumors and identify differentially expressed genes in cancer cells compared to normal cells.
## The Main Goal
In this project, I have selected a high-level article about employing scRNA-seq and bulk RNA-seq to observe the response to Glucocorticoids in human breast cancer cells. Then, I attempted to analyze the raw data (deposited in the [GEO database](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141834)) and rewrite the code exploiting R programming in order to generate plots and figures of the article.
<br> The main goal of this project is to examine how much I can analyze the sequencing data of cancer cells and derive valuable outcomes that can be employed for battling cancer.
## Brief Summary of The Article
First of all, this research article is a report on a high-level project, conducted by Hoffman et al. This article was published in the Communications Biology journal and the title of the paper is "single-cell RNA sequencing reveals a heterogeneous response to Glucocorticoids in breast cancer cells".
<br> Cellular heterogeneity is one of the main obstacles in the way of profiling and treatment of cancers. The main goal that Hoffman and his colleagues pursue in this study is to investigate the cellular heterogeneity in the transcriptional response of human breast cancer cells to glucocorticoids. To tackle this, they opt to detect transcriptional response, employing scRNA-seq across 18 hours of dexamethasone (Dex) treatment. Ultimately, scRNA-seq revealed significant cellular heterogeneity in response to Dex treatment from the glucocorticoid receptor target genes. Furthermore, the results unfold that there is remarkable variability between individual cells in the transcriptional response to hormones.
## Softwares and Packages
R programming (version 4.1.2) has been used for statistical computing and plotting figures. In addition, the Seurat R package (version 3.0.2) has been employed to perform normalizing, scaling, principal component analysis, and differential expression analysis on the single-cell data. Furthermore, the ggplot2 R package (version 3.3.5) has been exploited to plot the results of the Seurat package. Additionally, Limma (version 3.50.3) and edgeR (version 3.36.0) R packages have been used for finding differentially expressed genes in the bulk RNA-seq data. 
## Data Analysis and Plotting Figures
### 1- Basic Data Preparation
In the first step, the raw data is called and converted it into a Seurat object. 
```
dex <- read.delim("Data/GSE141834_scRNAseq_rawCounts.txt", header = T, row.names = 1, sep = "\t")
dex <- CreateSeuratObject(dex,project = "DEX")
```
Then, the data is normalized using the NormalizeData function and the cells with the percentage of mitochondrial genes higher than 5% is removed by the following line of codes:
```
dex$percent.mt <- PercentageFeatureSet(dex, pattern = "^MT-")
dex <- subset(dex, subset = percent.mt < 5)
```
After that, we can plot the some of the figures of the article by means of ggplot2 package such as the following graphs:
<img src="/Plots/Fig 1f.png" alt="Figure 1f" class="center" width="250">
<img src="/Plots/Fig 1g.png" alt="Figure 2g" class="center" width="250">
<img src="/Plots/Fig 1h.png" alt="Figure 2h" class="center" width="250">
### 2- Cell Cycle Scoring
In the next phase, the cell cycle scoring is performed employing the cyclone function of the scran Package. However, in the first place, the gene names have to be converted into standard ones in order that has the best scoring. After finding the new gene names, the scoring is performed by the following lines of code:
```
dex.sce <- as.SingleCellExperiment(dex)
cycle.scores <- scran::cyclone(dex.sce, hs.pairs, gene.names = new.rownames$gene_name)
```
### 3- Scaling and Principal Component Analysis (PCA)
In this step, the top 500 variable genes is scaled and centered exploiting negative binomial model. The impacts of transcript counts, mitochondrial percent, and cell-cycle scores has been regressed-out.
```
dex <- FindVariableFeatures(dex, nfeatures = 500)
dex <- ScaleData(dex, model.use = "negbinom", vars.to.regress = c("percent.mt", "G1","S","G2M", "nCount_RNA"))
```
Then, PCA is performed on the data using the 30 principal components and the results is plotted on a scatter plot:
```
dex <- RunPCA(dex, verbose = F, npcs = 30, seed.use = 119)
DimPlot(dex, reduction = "pca", group.by = "Hours_Dex", pt.size = 1.5) +
  scale_y_continuous(breaks = seq(-15,10,5))+
  theme(legend.title = element_text(size = 20),
        axis.title = element_text(size = 17),
        axis.text = element_text(size = 17),
        legend.text = element_text(size = 17)) + 
  guides(color= guide_legend(title = "Hrs Dex"))
```
<img src="/Plots/Fig 3a.png" alt="Figure 3a" width="400">
Likewise, tSNE and UMAP are performed on the data using the top 16 principal components.
### 4- Clustering the Data
The next phase comprises of finding the nearest neighbours using the top 16 principal components. After that, the data is clustered into 7 clusters using a proper resolution.
```
dex <- FindNeighbors(dex, reduction = "pca", dims = 1:16)
dex <- FindClusters(dex, resolution = 0.45)
```
Then, the number of cells in each cluster in different Dex treatment timepoints is plotted in the following histogram:
<img src="/Plots/Fig 3g.png" alt="Figure 3g" width="700">
### 5- Differentially Expressed Genes (DEGs)
After clustering the data, the differentially expressed genes is found by comparing each dex-treatment timepoint with the untreated cells. 
For this matter, we exploit the MAST test with a fold change cutoff of 1.25, an adjusted p-value of 0.01, and excluding genes detected in fewer than 10% of cells.
Ultimately, the following number of DEGs is detected in each Dex trwatment timepoint:
<img src="/Plots/Supplemental Fig 2b.JPG" alt="Supplemental Figure 2b" width="400">
Likewise, the DEGs of bulk RNA-seq is found by means of Limma-voom package. DEGs are called and filtered for each Dex-treatment timepoint with a fold change cutoff of 1.5 and adjusted p-value of 0.05.
The following number of DEGs is detected in each Dex trwatment timepoint:
<img src="/Plots/Supplemental Fig 2a.JPG" alt="Supplemental Figure 2a" width="400">
### 6- Ratio of Responding Genes (RRG)
Finally, the number Dex target genes showed a response in each cell is determined to calculate the ratio of responding genes.
Firstly, the mean log-scaled expression level and standard deviation is calculated for each DEG in untreated cells (using only non-zero values). A gene is "responsive" if it was expressed greater than one SD above the mean untreated level or more than one SD below the mean for downregulated genes.
<img src="/Plots/Fig 4d.png" alt="Figure 4d" width="400">
