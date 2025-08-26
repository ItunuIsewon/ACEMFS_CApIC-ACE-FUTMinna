# Introduction to R

### Biostatistics and Data Interpretation in R for Mycotoxicology
### Author: Dr Itunu I.M

### 📥 **Dataset:** Download the file [here](https://drive.google.com/file/d/1umPgNdUgH6e023l441B4fFO-MNk29ZEt/view?usp=drive_link)

Download the second file [here](https://drive.google.com/file/d/1vLykOWzcUNRORkfmQvFEXOI3-ouWMaXz/view?usp=drive_link)

### Introduction

Welcome to the workshop on **Biostatistics and Data Interpretation in R** with applications in **Mycotoxicology**. This training will guide participants from the basics of R programming to the application of statistical techniques and visualization methods in analyzing biological data related to mycotoxins.

---

# Getting Started with R and RStudio

```{r}
# Print Hello World
print("Hello, FUT Minna! Welcome to Biostatistics with R")
```

```{r}
# Basic arithmetic operations in R
2 + 3
10 / 2
sqrt(16)
```

```{r}
# Assigning variables
x <- 10
y <- 5
z <- x + y
z
```
In R, a variable is essentially a name that stores a value (data). Think of it like a container or a label that points to something you want to use later, such as a number, text, vector, dataframe, or even a function.

```{r}
# Creating vectors
toxins <- c("Aflatoxin", "Ochratoxin", "Fumonisin", "Zearalenone")
levels <- c(23, 45, 12, 30)
toxins
levels
```
In R, a vector is the most basic data structure — it is simply a sequence of elements of the same type (all numbers, all characters, all logical values, etc.).

---

### Loading and Cleaning Biological Data in R

```{r}
# Load necessary libraries
library(readr)
library(dplyr)

# Import CSV file of mycotoxin levels
myco <- read_csv("mycotoxin_data.csv")
head(myco)

```

```{r}
# Checking summary statistics
summary(myco)
```

```{r}
# Check missing values
colSums(is.na(myco))
```

```{r}
# Remove missing values
myco_clean <- na.omit(myco)

```

```{r}
# Remove missing values
colSums(is.na(myco_clean))

```

```{r}

# Rename columns for consistency
myco_clean <- myco_clean %>%
  rename(
     Aflatoxin_B1= "AFB1",
    Fumonisin_B1  = "FB1"
  )
```

---

### Basic Statistical Tests

### T-test: Comparing Aflatoxin levels in two food groups

```{r}
maize <- subset(myco_clean, Crop == "Maize")$AFB1
rice  <- subset(myco_clean, Crop == "Rice")$AFB1

t.test(maize, rice)

```
🔎 **Interpretation of Results**

Groups compared

maize vs rice (likely AFB1 or FB1 concentrations in your dataset).

Test statistic

t = 0.3124 → the difference between the two group means, scaled by variability, is very small.

Degrees of freedom (df)

df = 37.866 → Welch’s test adjusts df since variances might not be equal.

p-value

p-value = 0.7564 → This is much greater than 0.05, so there is no statistically significant difference between maize and rice toxin levels in this sample.

Confidence interval

95% CI: -4.88 to 6.66 → Since this interval includes 0, the difference could be negative or positive. Again, no evidence of a true difference.

Means

mean of x (maize) = 15.60

mean of y (rice) = 14.71
→ The averages are close, consistent with the non-significant result.

### ANOVA: Testing differences across multiple mycotoxin groups

```{r}
aflatoxin <- c(12, 15, 14, 18, 20, 22)
ochratoxin <- c(8, 10, 9, 11, 13, 12)
fumonisin <- c(25, 27, 26, 30, 28, 29)

myco_levels <- data.frame(
  value = c(aflatoxin, ochratoxin, fumonisin),
  group = rep(c("Aflatoxin", "Ochratoxin", "Fumonisin"), each = 6)
)

anova_result <- aov(value ~ group, data = myco_levels)
summary(anova_result)
```
group

This is your factor (e.g., crop type).

Df = 2 → There are 3 groups (because df = number of groups − 1).

Sum Sq = 885.8 → Variation explained by group differences.

Mean Sq = 442.9 → Average variation per group.

F value = 61.61 → Very large F statistic → strong evidence of differences.

Pr(>F) = 5.84e-08 → This is essentially 0.0000000584, way below 0.05.
✅ **Conclusion**: Highly significant difference between at least one pair of groups.

Residuals

Df = 15 → Total observations − number of groups.

Sum Sq = 107.8 → Variation not explained by the group (within-group variance).

Mean Sq = 7.2 → Average within-group variance.

Significance codes

*** means p < 0.001 → extremely strong evidence against the null hypothesis.

```{r}
#Question: Are there regional differences in Fumonisin B1 levels?
anova_model <- aov(FB1 ~ Region, data = myco_clean)
summary(anova_model)

```

---

### Data Visualization in Mycotoxicology

### Boxplots of Mycotoxin Levels
📊 What is a Box Plot?

A box plot (or “box-and-whisker plot”) is a way to summarize the distribution of a dataset.
It shows:

Median (Q2) → the central line in the box.

Interquartile Range (IQR = Q3 − Q1) → the box itself, showing the middle 50% of data.

Whiskers → extend to the smallest and largest values within 1.5×IQR.

Outliers → points beyond the whiskers.

🔎 When & Why to Use Box Plots

To compare distributions across groups (e.g., crops, regions, toxins).

To spot outliers and extreme contamination values.

To summarize variability (spread, skewness).
```{r}
myco_data <- read.csv("mycotoxin_data2.csv", header = TRUE)

boxplot(myco_data$Aflatoxin, myco_data$Ochratoxin, myco_data$Fumonisin,
        names = c("Aflatoxin", "Ochratoxin", "Fumonisin"),
        col = c("red", "blue", "green"),
        main = "Distribution of Mycotoxin Levels in Samples",
        ylab = "Concentration (ppb)")
```

### Histograms
📊 What is a Histogram?

A histogram is a type of plot that shows how data are distributed across intervals (called bins).

The x-axis represents ranges of values (intervals).

The y-axis represents the frequency (or proportion) of data points in each interval.

Unlike a bar chart (which compares categories), a histogram is specifically for continuous (numeric) data.

Think of it as a way to see the shape of your data distribution.

🔎 When do we use Histograms?

You use histograms when you want to:

Understand the distribution of a continuous variable

Is it normally distributed? Skewed left or right? Multimodal?

Check for outliers or extreme values

A long tail may suggest unusual values worth investigating.

Assess spread and variability

Do values cluster tightly or spread widely?

Compare distributions between groups (using faceted histograms or overlayed densities).

Decide on statistical methods

Many tests assume normal distribution; histograms help check assumptions.

✅ **Why use Histograms?**

**Quick insight:** They let you see if most values fall in a certain range.

**Exploratory data analysis (EDA):** Often the very first step when examining new data.

**Quality control:** Useful in fields like manufacturing or labs to monitor measurements.

**Risk assessment:** In toxicology/mycotoxicology, histograms help reveal if contamination is concentrated in a few samples or spread widely.

🌍** Context Examples - Mycotoxicology**

Plot a histogram of AFB1 concentrations across all maize samples.

You might see that most samples are <20 µg/kg, but a few exceed 100 µg/kg (long right tail = skewed distribution).

```{r}
hist(myco_data$Aflatoxin, col = "red", main = "Aflatoxin Distribution",
     xlab = "Concentration (ppb)", breaks = 5)
```
```{r}
library(ggplot2)

ggplot(myco_clean, aes(x = AFB1, fill = Crop)) +
  geom_histogram(bins = 5, alpha = 0.6) +
  labs(title = "Distribution of Aflatoxin B1 Levels by Crop",
       x = "Aflatoxin B1 (µg/kg)", y = "Count")

```

### Density Plot
📊 What is a Density Plot?

A density plot is like a smoothed version of a histogram.

Instead of showing counts per bin, it estimates the probability density of a continuous variable.

The area under the curve = 1 (100% of the data).

It helps reveal the shape of the distribution (e.g., bell-shaped, skewed, multimodal).

Think of it as a way to see the “outline” of the histogram without the jagged bins.

🔎 **When to Use a Density Plot**

When you want to focus on the distribution shape rather than exact counts.

When comparing multiple groups — overlaying density curves makes comparison easier than overlapping histograms.

For large datasets, density plots give a clearer picture than histograms.
```{r}
plot(density(myco_data$Ochratoxin), main = "Ochratoxin Density Plot", col = "blue", lwd = 2)
```

### Scatter Plot

```{r}
plot(myco_data$Aflatoxin, myco_data$Fumonisin,
     xlab = "Aflatoxin (ppb)", ylab = "Fumonisin (ppb)",
     main = "Scatter Plot: Aflatoxin vs Fumonisin",
     col = "purple", pch = 19)
```

### Correlation Heatmap
📊 What is a Correlation Heatmap?

A correlation heatmap is a visual tool used to show the strength and direction of relationships between multiple variables in a dataset.
It is based on the correlation coefficient (r), which ranges from −1 to +1:

+1 → perfect positive correlation (variables increase together).

0 → no correlation (no linear relationship).

−1 → perfect negative correlation (one increases, the other decreases).

The heatmap displays these correlations as a colored grid, where colors represent the magnitude and direction of correlation.

🔎 **When & Why to Use a Correlation Heatmap**

✅ To quickly identify relationships between variables (e.g., toxins measured in food samples).
✅ To detect multicollinearity before building statistical or machine learning models.
✅ To explore patterns and clusters of highly related variables.
✅ To guide feature selection and interpretation of experimental results.

```{r}
library(corrplot)
corr_matrix <- cor(myco_data[,2:4])
corrplot(corr_matrix, method = "color", addCoef.col = "white")
```
### Bar plot (mean per group)
📊 **What is a Bar Plot (Mean per Group)?**
A bar plot is a chart that displays the mean (average) of a variable across different groups. Each bar represents the group’s average value.

It shows:

Bars → the mean of each group.

Categories → groups being compared (e.g., crops, regions, treatments).

🔎 **When & Why to Use Bar Plots (Mean per Group)**

To compare average levels across groups (e.g., mean aflatoxin concentration in maize vs. groundnut).

To simplify data by focusing on central tendency rather than full distribution.

Useful when raw data is noisy but group differences in means matter.

Great for presentations where a clear group-level summary is needed.

⚠️** Note: **Unlike box plots, bar plots hide distribution details and outliers, so they are best when you only want to highlight group averages.
```{r}
mean_data <- myco_clean %>%
  group_by(Crop) %>% summarise(mean_AFB1 = mean(AFB1), .groups="drop")

ggplot(mean_data, aes(Crop, mean_AFB1, fill=Crop)) +
  geom_col() +
  geom_text(aes(label=round(mean_AFB1,1)), vjust=-0.5) +
  labs(title="Mean Aflatoxin B1 by Crop", y="AFB1 (µg/kg)")
```

### Line plot (trend over time)
📈 **What is a Line Plot?**
A line plot (or time series plot) is a way to visualize how a variable changes over time or across an ordered sequence.
It connects individual data points with lines, making trends easier to see.

It shows:

X-axis → usually time (days, months, years, seasons, storage duration).

Y-axis → values of the variable being measured (e.g., toxin concentration).

Lines → show upward or downward trends.

Multiple lines → allow comparison across groups (e.g., different toxins, regions).

🔎** When & Why to Use Line Plots**

To observe trends (increasing, decreasing, stable) over time.

To compare trajectories of multiple groups (e.g., Aflatoxin vs. Fumonisin levels).

To detect seasonal or temporal patterns in contamination.

To visualize the effect of interventions (e.g., storage method, treatment).

```{r}
time_data <- tibble(
  Month = 1:6,
  AFB1 = c(12, 20, 35, 50, 70, 90)
)

ggplot(time_data, aes(Month, AFB1)) +
  geom_line(linewidth=1) + geom_point(size=2) +
  labs(title="AFB1 Increase During Storage", y="AFB1 (µg/kg)")
```

### Violin plot (distribution + density)
🎻 **What is a Violin Plot?**
A violin plot combines features of a box plot and a density plot:

Like a box plot, it shows median, quartiles, and range.

Like a density plot, it shows the distribution shape (wider sections = more frequent values, thinner sections = fewer).

It looks like a “mirrored” density plot around a central line, resembling a violin 🎻.

🔎 **When & Why to Use Violin Plots**

To compare the distribution of a variable across groups.

To see spread + central tendency + density shape all in one plot.

Useful when data is not normally distributed (because the density shape reveals skewness/multi-modality).
```{r}
ggplot(myco_clean, aes(Crop, AFB1, fill=Crop)) +
  geom_violin(trim=FALSE, alpha=.7) +
  geom_boxplot(width=.12, fill="white") +
  labs(title="Distribution of AFB1 by Crop", y="AFB1 (µg/kg)")
```

```{r}
ggplot(myco_clean, aes(Region, AFB1, fill=Region)) +
  geom_violin(trim=FALSE, alpha=.7) +
  geom_boxplot(width=.12, fill="white") +
  labs(title="Distribution of AFB1 by Crop", y="AFB1 (µg/kg)")
```
### Density plot (smooth distribution)
📊 **What is a Density Plot?**
A density plot is a smooth curve that estimates the probability distribution of a dataset.
It is essentially a smoothed version of a histogram, using kernel density estimation (KDE).

It shows:

The shape of the distribution (e.g., unimodal, bimodal, skewed).

Peaks (modes) where data are concentrated.

The spread of the data across values.

🔎 **When & Why to Use Density Plots**

To visualize the underlying distribution of continuous data.

To compare distributions between groups (e.g., mycotoxin levels in maize vs. wheat).

To see whether data are normally distributed or skewed.

Better than histograms when you want a smooth and interpretable shape.

👉 Great for spotting whether contamination levels are centered around low values, or if multiple peaks suggest subgroups in the data.

```{r}
ggplot(myco_clean, aes(AFB1, fill=Crop)) +
  geom_density(alpha=.5) +
  labs(title="Density of AFB1 by Crop", x="AFB1 (µg/kg)")
```

### Stacked bar (proportions by category)
📊 **What is a Stacked Bar Plot?**
A stacked bar plot shows the composition of categories within groups by stacking segments on top of each other.
Instead of showing only totals, it breaks them down into subgroups.

It shows:

Group totals → bar height (e.g., total samples per crop).

Subcategories → stacked colors (e.g., toxin types within each crop).

Relative proportions → can be displayed as raw counts or percentages.

🔎 **When & Why to Use Stacked Bar Plots**
✔ To visualize composition of groups (e.g., toxin types per crop).
✔ To compare proportions across categories.
✔ To see distribution trends across groups.
✔ Great for categorical breakdown of survey, contamination, or frequency data.
```{r}
cats <- myco_clean %>%
  mutate(Category = cut(AFB1, breaks=c(-Inf,20,50,Inf),
                        labels=c("Low (<20)","Medium (20–50)","High (>50)"))) %>%
  count(Crop, Category)

ggplot(cats, aes(Crop, n, fill=Category)) +
  geom_col() +
  labs(title="AFB1 Risk Categories by Crop", y="Sample count")+
  scale_y_continuous(labels = scales::percent)
```



### Boxplot + jitter (summary + raw points)
📦 **What is a Boxplot + Jitter?**
A boxplot summarizes the distribution of a dataset:

The box shows the interquartile range (IQR = 25th–75th percentile).

The line inside = median.

Whiskers extend to min/max within 1.5 × IQR.

Dots outside = outliers.

Jittered points overlay the raw data, showing the spread of individual observations.

This helps avoid overplotting by spreading out points horizontally.

🔎 **When & Why to Use**

To combine summary statistics (boxplot) with individual-level data (jitter).

Useful when you want to:

Show group differences (e.g., toxin levels across crops).

Detect outliers and variability.

Avoid misleading conclusions from boxplots alone.
```{r}
ggplot(myco_clean, aes(Region, FB1, fill=Region)) +
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  geom_jitter(width=.15, alpha=.6) +
  labs(title="FB1 by Region (with samples)", y="FB1 (µg/kg)")
```

### 4.J Regional comparison (simple “geo” bar)

```{r}
region_mean <- myco_clean %>%
  group_by(Region) %>% summarise(mean_AFB1 = mean(AFB1), .groups="drop")

ggplot(region_mean, aes(Region, mean_AFB1, fill=Region)) +
  geom_col() +
  labs(title="Regional Mean AFB1", y="AFB1 (µg/kg)")
```


### Principal Component Analysis (PCA)
📊 **What is PCA (Principal Component Analysis)?**
PCA is a dimensionality reduction technique that transforms a large set of correlated variables into a smaller set of uncorrelated components while preserving as much variation in the data as possible.

Principal Components (PCs): New variables formed as linear combinations of the original variables.

PC1: Captures the largest variance in the dataset.

PC2: Captures the second-largest variance, orthogonal to PC1.

Usually, the first 2–3 PCs explain most of the structure in the data.

🔎 **When & Why to Use PCA?**

To reduce dimensionality while keeping most of the information.

To visualize high-dimensional data in 2D or 3D plots.

To identify patterns and clusters in complex datasets (e.g., samples, crops, toxins).

To remove redundancy when variables are highly correlated.
```{r}
myco_scaled <- scale(myco_data[,2:4])
pca_result <- prcomp(myco_scaled, center = TRUE, scale. = TRUE)
summary(pca_result)

biplot(pca_result, col = c("red", "blue"))
# Scree plot from base R
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
plot(var_explained, type = "b",
     xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     main = "Scree Plot")

# Cumulative variance
lines(cumsum(var_explained), type = "b", col = "red")

```

---

### Conclusion

This workshop introduced participants to:

- Basics of R programming  
- Data handling and cleaning in R  
- Statistical tests applied to biological/mycotoxin data  
- Visualization techniques (boxplots, histograms, scatterplots, density plots, correlation heatmaps, PCA)  

These tools provide a foundation for analyzing complex datasets in **mycotoxicology research**, such as understanding contamination levels, testing for differences across food groups, and exploring multivariate data structures.


