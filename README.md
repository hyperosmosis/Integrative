# Integrative Analysis for Microbiome Data Analysis

Metabolic dysregulation and alterations have been linked to various diseases and conditions. Innovations in high-throughput technology now allow rapid profiling of the metabolome and metagenome — often the gene content of bacterial populations -– for characterizing metabolism. Due to the small sample sizes and high dimensionality of the data, pathway analysis (wherein the effect of multiple genes or metabolites on an outcome is cumulatively assessed) of metabolomic data is commonly conducted and also represents a standard for metagenomic analysis. However, how to integrate both data types remains unclear. Recognizing that a metabolic pathway can be complementarily characterized by both metagenomics and metabolomics, we propose a weighted variance components framework to test if the joint effect of genes and metabolites in a biological pathway is associated with outcomes.

Description of the files are below
* 'WK_Functions.R' includes all the helper functions for the method and requires installation of the 'SKAT' package.
* 'WK_Type1.R' includes code for running simulations on the assessing Type 1 error of the method
* 'WK_Analysis.R' includes code for analyzing and visualizing the results of the simulations
