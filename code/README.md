# **ADNI QPP Template Analysis**
This repository contains MATLAB scripts for analyzing **Quasi-Periodic Patterns (QPPs)** from the **ADNI dataset**. The analysis focuses on **waveform visualization, correlation analysis, and projection onto fMRI time series**.

## **Overview**
This repository includes three primary MATLAB scripts:
- **QPP Waveform Visualization:** Plots QPP waveforms across Neuromark networks.
- **QPP Correlation Analysis:** Computes correlation matrices across different QPP templates and evaluates significance.
- **QPP Projection:** Projects QPP templates onto fMRI time series and evaluates occurrences over time.

Below are example outputs from each section:

---

## **1️⃣ QPP Waveform Visualization**
📜 **Script:** `m1_QPP_waveform_visualization_neuromark.m`  
📊 **Purpose:** Visualize QPP waveforms across different brain networks and subdomains.

### **Example Outputs**
- **QPP Waveforms Across Networks**  
  *Visualization of waveform patterns across multiple Neuromark networks.*  
  ![QPP Waveforms](plots/m1_qpp_waveforms.png)

- **Network-Specific Coloring of QPPs**  
  *Each network is color-coded to show activity over time.*  
  ![QPP Network Coloring](plots/m1_qpp_network_coloring.png)

- **HC vs TN Correlation**  
  *Comparison of Quasi-Periodic Patterns between the HC and TN networks.*  
  ![HC-TN Correlation](plots/m1_hc_tn_correlation.png)

- **Network Correlation Matrix**  
  *Shows correlations between different Neuromark networks based on QPP signals.*  
  ![Network Correlation Matrix](plots/m1_network_correlation_matrix.png)

- **Subdomain Correlation Matrix**  
  *A correlation matrix focused on subdomains of the Neuromark networks.*  
  ![Subdomain Correlation Matrix](plots/m1_subdomain_correlation_matrix.png)

---

## **2️⃣ QPP Correlation Analysis**
📜 **Script:** `m2_QPP_corr_analysis.m`  
📊 **Purpose:** Compute **correlation matrices** between different QPP templates and evaluate network significance.

### **Example Outputs**
- **QPP Template Correlation Matrix**  
  *Compares correlation across all QPP templates.*  
  ![Template Correlation Matrix](plots/m2_qpp_template_correlation_matrix.png)

- **QPP Template Correlation Analysis**  
  *Computation of self-correlation differences between QPP templates.*  
  ![Template Correlation Analysis](plots/m2_qpp_template_correlation_analysis.png)

- **QPP Correlation Comparison**  
  *Comparing correlation matrices between different QPP templates.*  
  ![QPP Correlation Comparison](plots/m2_qpp_correlation_comparison.png)

- **QPP Difference Correlation Matrix**  
  *Visualizes the difference between QPP correlation matrices.*  
  ![QPP Difference Correlation](plots/m2_qpp_difference_correlation_comparison.png)

- **QPP Full Correlation Matrix Comparison**  
  *Comparison of full correlation matrices across templates.*  
  ![QPP Full Correlation](plots/m2_qpp_full_correlation_comparison.png)

- **ICN Significance in QPP Correlation**  
  *Significance testing of Independent Component Networks (ICNs).*  
  ![ICN Significance](plots/m2_qpp_correlation_ICN_significance.png)

- **Network Variance Comparison**  
  *Comparison of network variance across templates.*  
  ![Network Variance](plots/m2_network_variance_comparison.png)

- **Boxplot of Mean Correlation Values**  
  *Boxplot visualization for comparison of QPP template correlation means.*  
  ![Boxplot Means Visualization](plots/m2_boxplot_means_visualization.png)

---

## **3️⃣ QPP Projection**
📜 **Script:** `m3_QPP_projection.m`  
📊 **Purpose:** Project QPP templates onto **fMRI time series** and analyze **dwell times**.

### **Example Outputs**
- **Correlation Plot Over Time**  
  *Sliding window correlation of QPP projection onto fMRI data.*  
  ![Correlation Plot](plots/m3_correlation_plot.png)

- **Correlation Histogram**  
  *Histogram of correlation values across sliding windows.*  
  ![Correlation Histogram](plots/m3_correlation_histogram.png)

- **Correlation with Dwell Time Markers**  
  *Adding dwell time markers to the correlation plot.*  
  ![Correlation with Dwell](plots/m3_correlation_plot_with_dwell.png)

- **Dwell Time Histogram**  
  *Histogram of time intervals between significant correlation events.*  
  ![Dwell Time Histogram](plots/m3_dwell_time_histogram.png)

---

## **4️⃣ Python Data Visualization**
📜 **Script:** `python_data_visualization.ipynb`  
📊 **Purpose:** This **standalone Jupyter Notebook** provides additional visualization for the **longitudinal dataset**.  
It focuses on:
- **Graphing projection values** over time.
- **Visualizing changes** in QPP occurrences across different conditions.
- **Comparing QPP activity** across different ADNI cohorts.

### **Example Outputs**
- **Longitudinal Transition Point Graph**  
  *Visual representation of longitudinal changes in QPP projection values.*  
  ![Longitudinal Transition Point Graph](../figures/figure_2_v5.PNG)

- **QPP Projection Results**  
  *QPP projection analysis across different timepoints and conditions.*  
  ![QPP Projection Results](../figures/figure_6_v3.PNG)
