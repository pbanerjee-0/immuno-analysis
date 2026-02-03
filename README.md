# Computational Analysis of Immune Cell Dynamics

## Project Overview
This repository contains a specialized computational pipeline developed during my tenure at the Starzl Institute (UPMC) to analyze high-resolution longitudinal imaging data of T-cells (CD8+, CD4+, and Tregs). The primary objective is the quantification of calcium signaling dynamics—a critical indicator of immune response activation.

## The Challenge: Dynamic Baseline Normalization
In biological datasets, calcium spike intensity and average concentrations fluctuate significantly across different cell tracks and experiments. To ensure accurate analysis, the definition of a "calcium spike" must be adaptive rather than static. This pipeline implements a statistical approach to establish dynamic thresholds for every individual dataset.

## Technical Workflow & Methodology

### 1. Signal Pre-processing & Thresholding
* **Linear Regression Modeling:** Implemented Python-based regression scripts to fit raw fluorescence intensity data to a curve.
* **Dynamic Baseline Establishment:** The fitted curve is used to calculate a specific calcium spike threshold for each individual "movie" or dataset, accounting for background noise and varying intensity baselines.

### 2. Feature Extraction & Peak Detection
* **Spike Quantification:** Using the derived thresholds as parameters, the pipeline automates the detection and counting of calcium spikes across thousands of timepoints.
* **Kinetic Analysis:** Leveraged `SciPy` peak detection algorithms to determine the frequency and magnitude of signaling events.
* **Area Under the Curve (AUC):** Utilized the trapezoidal rule to calculate the duration and total integrated intensity of each calcium spike, providing a measure of total signaling "dosage."

### 3. Biological Context & Output
* **Multivariate Analysis:** The pipeline correlates calcium signaling frequency with APC (Antigen-Presenting Cell) contact duration.
* **Data Export:** Coefficients and statistical summaries are exported for high-level visualization in GraphPad Prism and Matplotlib to analyze the effects of mTOR inhibitors on immune cell dynamics.

## Tech Stack
* **Language:** Python
* **Data Manipulation:** Pandas, NumPy
* **Scientific Computing:** SciPy (Signal Processing)
* **Visualization:** Matplotlib
* **Statistical Analysis:** Linear Regression

## Impact
This pipeline standardizes the analysis of immune cell dynamics, allowing for the rapid processing of massive datasets (3,000–10,000 timepoints per cell) while maintaining high sensitivity to biological variation.
