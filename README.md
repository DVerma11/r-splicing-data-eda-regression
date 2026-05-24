# Exploratory Data Analysis and Predictive Modeling with Synthetic Splicing Data

This R project performs exploratory data analysis, visualization, correlation analysis, and predictive modeling using a synthetic gene-splicing dataset.

## Objective
The goal is to examine relationships between three splicing factors and a splicing event, then build a linear regression model to predict the splicing event.

## Methods
- Summary statistics
- Missing value check
- Outlier detection using boxplots
- Correlation matrix using Pearson correlation
- Histograms with normal curves
- Scatterplots with regression lines
- Linear regression modeling
- Model evaluation using R-squared and Mean Squared Error

## Key Findings
- No missing values were found.
- Two outliers were detected.
- SplicingFactor1 was positively associated with the splicing event.
- SplicingFactor2 was negatively associated with the splicing event.
- SplicingFactor3 showed little to no relationship.
- Linear regression showed moderate to strong predictive performance.

## Requirements
```r
install.packages(c("ggplot2", "Hmisc", "pastecs"))
