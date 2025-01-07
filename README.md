# Codon Adaptiveness Index (CAI) calculation
Decoding Expression Potential - A Python Program for Codon Adaptiveness Index (CAI) Calculation

The Codon Adaptiveness Index (CAI) serves as a pivotal metric in understanding gene expression. This project
addresses the analysis of Codon Usage Bias in three strains of the Zika virus: PRVABC59,
ZIKV/Hu/Chiba/S36/2016 (ChibaS36), and ZIKV/Hu/NIID123/2016 (NIID123). Codon usage bias can provide
insights into viral evolution, adaptation, and host interactions. This is a Python program that calculates
Relative Synonymous Codon Usage (RSCU), Codon Adaptation Index (CAI), and facilitates comparison across
samples through visualizations such as bar graphs and heatmaps. 

The program takes three gene samples as input and efficiently compares their CAI values, 
pinpointing the sample with the higher expression potential.
Leveraging the relative adaptiveness parameter or relative synonymous codon usage (RCSU), the program
computes the frequency of each codon relative to the most frequently used synonymous codon within a set of
highly expressed genes. Subsequently, the CAI for each coding sequence is derived as the geometric mean of
their respective RCSU values. 

This program offers a readily accessible and efficient tool for researchers to
rapidly predict highly expressed genes, compare expression potential across samples and easily integrate the
program into existing bioinformatics pipelines. This project contributes to advancing gene expression prediction,
aiding research across various fields, including synthetic biology, genetic engineering, biotechnology and
beyond.
