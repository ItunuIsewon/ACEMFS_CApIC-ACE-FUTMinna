## GROUP 1 PROJECT

**Author:** Itunuoluwa Isewon PhD     

**Email:** itunu.isewon@covenantuniversity.edu.ng

## Phylogenetic Analysis of Fungal 28S Sequences

### Background

The 28S region of rDNA is widely used as a DNA barcode for fungi. In this project, students will analyse 28S sequences from five fungal species, along with a set of unknown fungal sequences. The aim is to use phylogenetic approaches to determine which species the unknown belongs to, and to compare different tree-building methods.

### Objectives

•	Perform multiple sequence alignment (MSA) of fungal 28S sequences. 
•	Construct phylogenetic trees using Neighbour Joining (NJ), Maximum Likelihood (ML), and Minimum Evolution (ME) methods.
•	Compare phylogenetic tree topologies across methods.
•	Identify the most likely species for the unknown sequences.

### Dataset Provided

### 📥 **Dataset:**

FASTA file containing: 

•	28S sequences from 5 known fungal species
•	Several unknown 28S sequences to classify

Download dataset [here](https://drive.google.com/file/d/1Du-DAWhij6LZp5tCMSdW5byBG34W07t2/view?usp=sharing)

### Tasks

1.	Load the sequences Import the FASTA file into MEGA
2.	Multiple Sequence Alignment (MSA) Using ClustalW to align the 28S sequences.
3.	Save the alignment file. 
4.	Phylogenetic Tree Construction: Build trees using three methods;
•   Neighbor Joining (NJ)
•	 Maximum Likelihood (ML)
•   Minimum Evolution (ME)
6.	Visualize each tree and annotate species names.
   
### Analysis

•	Compare the positions of the unknown sequences in each tree. 
•	Check if all methods agree on the classification of unknowns. 
•	Note any differences in tree topology between methods.

### Expected Output
1.	A ClustalW alignment file
2.	Three phylogenetic trees (NJ, ML, ME)
3.	A presentation including: 
•   Figures of the phylogenetic trees 
•	 Identification of unknown sequences
•	 Discussion of agreement/disagreement among methods

## GROUP 2 PROJECT

## Transcriptomics (RNA-Seq Differential Expression Analysis)

### Introduction

In this mini-project, participants will perform differential gene expression analysis using RNA-seq data from the NCBI Gene Expression Omnibus (GEO) dataset GSE292521. Raw counts will be provided in a CSV format along with sample group information. This exercise will allow participants to gain hands-on experience with transcriptomics data analysis.

### Summary of GSE292521 Experiment

**Title:** Genomics and Transcriptomics of 3ANX (NX-2) and NX (NX-3) producing isolates of Fusarium graminearum

**Organism:** Fusarium graminearum

**Scope:** Analysis of 20 fungal isolates from different regions in Manitoba, characterizing both genomic variations and gene expression profiles linked to mycotoxin chemotypes (3ANX and NX). The data illuminate differential expression patterns related to pathogenicity and suggest the 3ANX chemotype may be more widespread in Canada than previously recognised.

**Platform:** Illumina NovaSeq 6000 sequencing.

### Objectives

By the end of this mini-project, participants will be able to:
•	Import RNA-seq count data and metadata into R.
•	Conduct differential expression analysis using DESeq2.
•	Visualize results using MA plots, heatmaps, and volcano plots.
•	Interpret biological significance of differentially expressed genes.

### 📥 **Dataset:**

**Dataset:** GSE292521 (Available on GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292521)
Provided Files:
+ Raw counts in CSV format
+ Metadata (sample groupings, e.g., control vs treatment)
  
Download dataset [here](https://drive.google.com/file/d/1icHEWAFlmqGRfxUaq6zIF9Sv9LUBNq3u/view?usp=sharing)

### Tasks

Participants are expected to carry out the following steps:
1.	Load count data and metadata into R.
2.	Normalize data using DESeq2.
3.	Conduct differential expression analysis.
4.	Generate visualization plots (volcano plot, heatmap).
5.	Identify top differentially expressed genes and their biological roles.
   
**Required R Packages**

The following R packages will be needed for this project:
- DESeq2 (for differential expression analysis)
- pheatmap (for heatmaps)
- EnhancedVolcano (for volcano plots)
- ggplot2 (for visualization)
- dplyr (for data manipulation)
  
### Expected Output

By completing this project, participants should produce:
•	A heatmap of top differentially expressed genes.
•	A volcano plot highlighting significant genes.
•	A list of significantly upregulated and downregulated genes.
•	A Group presentation showing the rationale, all plots, the top differentially expressed genes, their associated enriched KEGG pathways and Link findings to mycology/mycotoxicology relevance (e.g., stress response, secondary metabolite genes, toxin biosynthesis genes if present)

