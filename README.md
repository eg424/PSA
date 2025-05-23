# Probability, Statistics & Data Analysis
This repository contains MATLAB code used to complete the Probability and Statistics coursework. The coursework explores statistical testing, confidence intervals, normality testing, outlier detection, power analysis, and linear regression. It is divided into five numerical questions, each requiring analysis of experimental data generated in MATLAB.

## Tasks Summary

### Question 1 – One-Sample t-test
Assessed whether OD values using Arginine (A1) are significantly greater than a baseline of 0.35.
- Manual computation of test statistic and _p_-value.
- Comparison with MATLAB's built-in `ttest`.

### Question 2 – Two-Sample t-test with Unequal Variance
Compared OD means for Urea (U2) and Arginine (A2).
- Calculated confidence intervals for the difference in means.
- Used Welch's t-test for inference.

### Question 3 – Nonparametric Testing and Outlier Detection
Evaluated normality using Shapiro-Wilk tests for U3 and A3.
- Detected and removed outliers using MAD (Median Absolute Deviation).
- Applied Mann-Whitney U test (Wilcoxon rank-sum) post-cleaning.

### Question 4 – Paired Testing and Power Analysis
Analyzed control (C1, C2) and test (T1, T2) data across 350 students.
- Calculated proportions of significant/non-significant p-values.
- Estimated statistical power using `sampsizepwr`.
- Determined additional sample size needed for 90% power with α=0.01.

### Question 5 – Correlation and Linear Regression
Assessed the relationship between antibiotic concentration and OD.
- Calculated Pearson's correlation and regression parameters.
- Constructed confidence and prediction intervals.
- Estimated OD prediction interval for 2 µg/ml concentration.

## Figures

### Figure 1 – Boxplot for Urea vs Arginine
Displays spread and mean OD values with color-coded boxplots and overlays.

### Figure 2 – p-Values Distribution (log10 scale)
Scatter plot comparing -log10(p-values) for control and test data across students.

### Figure 3 – Regression with Confidence and Prediction Bands
Scatter plot with fitted regression line, confidence bands, and prediction intervals.

## Statistical Notes

- All tests used α = 0.05 unless specified otherwise.
- Where appropriate, both manual and MATLAB built-in methods were used.
- Outliers were identified using a 3×MAD rule.
- Figures exported at high resolution as per instructions.
