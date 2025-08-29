## Mini- Projects for the Introduction to Bioinformatics Workshop

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


## GROUP 3 PROJECT

## Mini Project: End-to-End Genomic Data Analysis of Fungal Isolates (Galaxy)

### Background
Whole-genome resequencing (WGS) enables variant discovery and comparative genomics across fungal isolates. In this project, students will process short-read FASTQ files from 5 known fungal species plus several unknown isolates. They will perform quality control, read trimming, reference-guided alignment, variant calling, and annotation in Galaxy.

### Objectives

●	Navigate public databases to find references, annotations, and reads (NCBI/ENA/SRA, Ensembl Fungi).
●	Execute a complete WGS pipeline in Galaxy: QC → trimming → alignment → post-processing → joint variant calling → filtering → annotation → core SNP matrix.

### Dataset Provided

A project folder containing:
•	FASTQ (paired-end) reads for:
1.	5 known fungal species (≥2 isolates each if available).
2.	Several unknown isolates to classify.

•	A metadata sheet (samples.tsv) with columns: sample_id, status (Known/Unknown), species (for known), SRA_accession (if applicable), library_layout, read_group, notes.
If you prefer, provide SRA accessions only; students will fetch reads inside Galaxy using “NCBI SRA Tools”.

### Tasks

1) Database Navigation (NCBI/ENA/SRA & Ensembl Fungi)
Goal: identify and download all inputs reproducibly.
 Steps (students document each):
1.	Reference genome:

○	Find a chromosome-level or best-available assembly for each known species on NCBI Assembly or Ensembl Fungi.

○	Download: reference.fasta and annotation.gff3 (or GTF).

○	Record assembly accession (e.g., GCA_XXXXXXXXX.X) and version.

2.	Reads (if not pre-provided):

○	Use SRA Run Selector to list runs, confirm paired-end, Illumina platform, and similar read length.

○	Note the SRA accessions for each sample and add to samples.tsv.

3.	Document database pages/screenshots + accessions in a short “Data Provenance” note.
2) Load Data into Galaxy
Goal: organise the project as reproducible Galaxy histories & collections.
 Steps:
●	Create a Galaxy History named Fungal_Genomics_Project_<YourName>.

●	Upload or fetch:

○	FASTQs (or use Get Data → NCBI SRA Tools: Fasterq-dump).

○	References (reference.fasta) + annotation.gff3 for each species.

●	Build Dataset Collections for paired reads (R1/R2).

●	Rename datasets with clear labels: SPC1_iso1_R1 / SPC1_iso1_R2, etc.
3) Quality Control (Galaxy)

**Tools & sequence (typical choices in parentheses):**
1.	FastQC on all raw reads.

2.	MultiQC to summarize FastQC results.

3.	Adapter/quality trimming (e.g., Trim Galore! or fastp):

○	Typical params: adapter auto-detect; quality cutoff 20; min length 50–70.

4.	FastQC (post-trim) → MultiQC (compare improvement).
 Output: MultiQC HTML reports for raw and trimmed reads.

4) Alignment & Post-Processing (Galaxy)
Reference choice:
●	Option A (simplest): use a single reference from the species you expect most unknowns to belong to.

●	Option B (rigorous): map each isolate to its species-specific reference, then combine variants in a species-aware manner. (Pick A for first run; B as bonus.)

**Tools & steps (BWA-MEM2 pipeline example):**
1.	BWA-MEM/MEM2: index reference.fasta; map paired reads → SAM.

2.	Samtools sort → BAM; Samtools index.

3.	Picard MarkDuplicates (or GATK MarkDuplicates).

4.	Alignment metrics:
Samtools flagstat and idxstats;
Optional: Qualimap BamQC;
bedtools genomecov for coverage summaries.
Output: deduplicated, indexed BAMs; metrics tables; coverage summaries.

5) Variant Calling & Joint Genotyping (Galaxy)
Two solid routes (pick one):
Route 1: bcftools mpileup/call
1.	bcftools mpileup (per sample) with -Ou -f reference.fasta.

2.	bcftools call (per sample) with -mv (variants only) → per-sample VCF.

3.	bcftools merge (multi-sample) to form a joint VCF.

Route 2: FreeBayes (population calling)
●	FreeBayes on a collection of BAMs to emit one multi-sample VCF.

Filtering (either route):
●	bcftools filter: depth (e.g., DP≥8–10), quality (QUAL≥30), genotype quality (GQ≥20), missingness (retain sites with ≥80% genotyped).

●	Optionally vcftools: --max-missing 0.8, --maf 0.01.

Output: a filtered multi-sample VCF.

6) Variant Annotation (Galaxy)
●	Build or select a database for your species:
SnpEff: create/select the fungal genome database (Ensembl Fungi GFF3).
Alternatively, VEP (if available) for functional consequence.

●	Run SnpEff on the filtered VCF → annotated VCF (*.snpeff.vcf) + summary HTML.

Output: functionally annotated VCF + summary.

7) Core SNP Matrix (Galaxy)
   
**Goal: derive a core SNP alignment.**

1.	Extract biallelic SNPs only (e.g., bcftools view -v snps -m2 -M2).

2.	Create alignment:
vcf2phylip (Galaxy wrapper) or SNP-sites to convert VCF → PHYLIP/FASTA SNP alignment.
Analysis
●	Data quality effects: Do low-coverage samples behave erratically? Reference bias?

●	Variant filters sensitivity: Show how stricter/looser DP/QUAL thresholds affect clustering.

●	Functional signals (optional): Are there species-specific HIGH-impact variants?

### Expected Output

1.	Galaxy artifacts
Two MultiQC reports (pre- and post-trim).
Deduplicated, indexed BAM files + alignment/coverage metrics.
Filtered multi-sample VCF and SnpEff-annotated VCF + HTML summary.
Core SNP alignment (FASTA/PHYLIP).

2.	Presentation (10–12 slides)
Figures: MultiQC screenshots, pipeline schematic.
Discussion of agreement/disagreement; sensitivity to filters; limitations.

3.	Reproducibility bundle
Galaxy Workflow export (.ga), History export (.tar), and a README with tool versions/parameters.

### Practical Tips & Parameter Hints
●	**Trim Galore!** : Q=20, min len=50–70; auto-adapters.

●	**Alignment:** BWA-MEM2 default; ensure read groups if needed.

●	**Filtering:** start with DP≥10, GQ≥20, QUAL≥30, max missing ≤20%; then explore sensitivity.

●	**SnpEff:** confirm genome build matches reference; document database/version.


