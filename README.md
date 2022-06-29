# Analysis of Single-Cell and Bulk RNA-seq Data by R
Rewriting the code of a paper about employing scRNA-seq and bulk RNA-seq for research on breast cancer. 
<br> This project has been done mainly by means of the Seurat package (version 3) in R.</br>
<br> The main challenge of this project was to rewrite the code without knowing the value of parameters that have been set in various functions. As a result, I had to run a function several times to find the best value for a parameter.
## Table of Contents
- [Introduction](#introduction)

## Introduction
Cancer is one of the leading causes of death in the world and it can affect people of all ages. Cancer cells are a result of genetic mutations caused by endogenous or environmental factors. The development of malignant tumors is a complex process that is mostly uncovered and it needs to make manifold efforts to unravel the main pathway that every normal cell passes until it turns into a cancer cell. Obviously, shedding light on this multi-stage process to decipher the hidden aspects of cancer cells' progression.
<br> Poor prognosis, the spread of cancer cells into distant organs due to metastasis, evasion from the immune system, and resistance to medications are among the major hurdles in the way of cancer therapy. For this reason, researchers in this field have to get familiar with cancer cells and tumor microenvironment at different stages in order to uncover the vulnerable points of the disease for designing brand-new medications. 
<br> Nowadays, manifold sequencing techniques and methods can make an important contribution to achieving this goal. Especially, single-cell RNA sequencing (scRNA-seq) and bulk RNA sequencing are featured methods that can help to investigate cellular heterogeneity of tumors and identify differentially expressed genes in cancer cells compared to normal cells.
## The Main Goal
<br> In this project, I have selected a high-level paper about employing scRNA-seq and bulk RNA-seq to observe the response to Glucocorticoids in human breast cancer cells. Then, I attempted to analyze the raw data (deposited in the [GEO database](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141834)) and rewrite the code exploiting R programming in order to generate plots and figures of the paper.
<br> The main goal of this project is to examine myself that how much I can analyze the sequencing data of cancer cells and derive 
