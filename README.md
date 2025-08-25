# Maximum Likelihood Data Analysis with ROOT

This repository contains a statistical data analysis project using the ROOT framework for analyzing radioactive decay time distributions. The project focuses on implementing various maximum likelihood estimation techniques and comparing binned and unbinned fitting methods.

## Project Structure

The project is organized into multiple parts:

- **Part 1**: Implementation of binned data analysis with Poisson likelihood fitting.
- **Part 2**: Implementation of unbinned maximum likelihood estimation and binned fits with various bin counts.
- **Part 3**: Monte Carlo simulations and chi-squared analysis.

## Key Features

- Implementation of analytical maximum likelihood estimators
- Comparison of binned vs. unbinned maximum likelihood methods
- Visualization of log-likelihood functions
- Error estimation using different statistical approaches
- Analysis of fitting performance with varying bin counts (3, 5, and 350 bins)
- Colorful and informative data visualizations

## Technologies Used

- C++ programming language
- ROOT data analysis framework (CERN)
- Statistical analysis techniques (Maximum Likelihood, Chi-squared, Poisson statistics)
- Data visualization

## Getting Started

To run these analyses, you need to have ROOT installed on your system. You can run the macros with:

```bash
root -l ML_Project_Part2.c
root -l ML_Project_Part3.C
```

## Output Examples

The repository includes several output plots:
- Log-likelihood function plots
- Binned fit histograms with 3, 5, and 350 bins
- Poisson-based fits
- Chi-squared vs. observation count analysis
- Monte Carlo lambda distribution
