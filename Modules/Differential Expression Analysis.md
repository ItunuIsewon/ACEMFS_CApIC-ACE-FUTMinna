# Differential Expression Analysis
## Author: Dr Itunu I.M

### 📥 **Dataset:** RNA-seq Analysis of GSE179477 with DESeq2
Download the file [here](https://drive.google.com/file/d/19BcWcMqddtN9YdMqGONyyHJgcvF8LUQw/view?usp=drive_link).

### Introduction

In this tutorial, we will perform a **differential expression analysis** on RNA-seq data from the dataset **GSE179477**.  
Download the file [here](https://drive.google.com/file/d/19BcWcMqddtN9YdMqGONyyHJgcvF8LUQw/view?usp=drive_link).
We will use the **DESeq2** package, which is widely used for RNA-seq analysis because it models count data using the negative binomial distribution.  

The goal is to demonstrate key steps:  
1. Loading data  
2. Quality check and normalization  
3. Running differential expression  
4. Visualizing results  

---

### Step 1: Load Libraries

We start by loading the libraries required for RNA-seq analysis.  

```{r setup, message=FALSE, warning=FALSE}
# Core packages
library(DESeq2)       # Differential expression analysis
library(ggplot2)      # Visualization
library(pheatmap)     # Heatmap plotting
library(dplyr)        # Data wrangling

```

### Step 2: Load the Data

We will use the read count matrix file: GSE179477_readcount.csv.

Each row = a gene

Each column = a sample

Values = raw read counts
```{r}
# Load data
counts <- read.csv("GSE179477_readcount.csv", row.names = 1)

# Preview the data
head(counts)
dim(counts)   # number of genes x samples

```
👉 **Explanation**: RNA-seq starts with raw counts. These must be normalized before comparisons across conditions.

### Step 3: Create Metadata (Sample Information)

DESeq2 needs sample metadata describing the experimental conditions.
You should create a dataframe that matches samples to conditions (e.g., control vs treatment).
```{r}
# Example metadata (replace with actual design of GSE179477)
sample_info <- data.frame(
  row.names = colnames(counts),
  condition = c("Control","Control","Control","Treatment","Treatment","Treatment")
)

sample_info

```
👉 **Explanation:** Each sample must have a label that tells DESeq2 which group it belongs to.

### Step 4: Create DESeq2 Dataset

We combine the count matrix and sample information into a DESeqDataSet.
```{r}
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = sample_info,
  design = ~ condition
)

dds
```
👉 **Explanation**: The design formula specifies the variables to test. Here we compare Treatment vs Control.

### Step 5: Pre-Filtering Low Count Genes

Genes with very low counts in all samples are not informative. We filter them out.

```{r}
dds <- dds[rowSums(counts(dds)) > 10, ]
dds
```
👉 **Explanation:** Removing low-count genes reduces noise and improves statistical power.

### Step 6: Run Differential Expression Analysis

We now run the DESeq2 pipeline.
```{r}
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("condition","Treatment","Control"))
head(res)
```
👉 **Explanation**: DESeq2 fits a negative binomial model and tests for differences between groups.
contrast specifies the comparison: Treatment vs Control.

Columns in DESeq2 Results

baseMean

The average normalized expression of a gene across all samples (both conditions/groups).

A higher value means the gene is overall more highly expressed in the dataset.

log2FoldChange (log₂FC)

The log₂ of the fold change between two groups (e.g., Treatment vs. Control).

Positive values → gene is upregulated in the treatment group.

Negative values → gene is downregulated in the treatment group.

Example: log₂FC = 1 → expression doubled; log₂FC = -1 → expression halved.

lfcSE (log fold change Standard Error)

The standard error of the estimated log₂ fold change.

Smaller values mean more confidence in the fold change estimate.

stat

The test statistic used for hypothesis testing (usually a Wald statistic in DESeq2).

Larger absolute values indicate stronger evidence for differential expression.

pvalue

The raw p-value from the statistical test.

Measures the probability that the observed difference happened by chance.

Not corrected for multiple testing.

padj (adjusted p-value / FDR)

The p-value adjusted for multiple testing using the Benjamini–Hochberg method (False Discovery Rate control).

This is the most reliable measure to decide significance.

A common cutoff: padj < 0.05 → significant differentially expressed gene.

**Example:**

| **Gene ** | **baseMean** | **log2FoldChange** |** lfcSE** | **stat ** | **pvalue ** | **padj**|
| ----- | -------- | -------------- | ----- | ----- | ------- | ------- |
| GeneA | 150.3    | 2.1            | 0.4   | 5.25  | 1.2e-07 | 3.5e-04 |
| GeneB | 75.8     | -1.3           | 0.6   | -2.16 | 0.031   | 0.12    |

GeneA: Highly expressed, significantly upregulated (log₂FC=2.1, padj < 0.05).

GeneB: Downregulated (log₂FC=-1.3), but not significant after correction (padj=0.12).

### Step 7: Data Visualization
MA Plot

Shows the relationship between gene expression (mean counts) and fold change.
```{r}
plotMA(res, ylim=c(-5,5))
```

Volcano Plot

Highlights significantly differentially expressed genes.
```{r}
res_df <- as.data.frame(res)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0, "Yes","No")

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
  geom_point(alpha=0.6) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  labs(title="Volcano Plot", x="Log2 Fold Change", y="-Log10 Adjusted p-value")

```

```{r volcano_plot, fig.width=10, fig.height=8}
# Assume results is your DESeq2 results dataframe
res_df <- as.data.frame(res)

# Add a column for significance category
res_df <- res_df %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  labs(title = "Volcano Plot of Differential Expression (DESeq2)",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())
```

```{r volcano_plot, fig.width=10, fig.height=8}
library(ggrepel)

# Convert DESeq2 results to dataframe
res_df <- as.data.frame(res)

# Add a column for significance category
res_df <- res_df %>%
  mutate(Significance = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "Upregulated",
    padj < 0.05 & log2FoldChange < 0 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 10 genes based on adjusted p-value
top_genes <- res_df %>%
  arrange(padj) %>%
  head(10)

# Volcano plot with annotations
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "grey")) +
  geom_text_repel(data = top_genes,
                  aes(label = rownames(top_genes)),
                  size = 4, box.padding = 0.4, point.padding = 0.3,
                  max.overlaps = Inf) +
  labs(title = "Volcano Plot of Differential Expression (DESeq2)",
       x = "log2(Fold Change)",
       y = "-log10(Adjusted p-value)") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank())

```


Heatmap of Top Genes
```{r fig.width=10, fig.height=8}
# Select top 30 significant genes
top_genes <- head(order(res$padj),30)

mat <- assay(vst(dds))[top_genes, ]
pheatmap(mat, annotation_col = sample_info, scale="row",
         main="Top 30 Differentially Expressed Genes")
```
👉 **Explanation:** Heatmaps show clustering of samples and expression patterns of significant genes.

### Step 8: Save Results

Finally, we save the results table.
```{r}
write.csv(as.data.frame(res), "GSE179477_DESeq2_results.csv")
```


