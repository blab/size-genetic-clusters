# Estimating the reproduction number and transmission heterogeneity from the size distribution of clusters of identical pathogen sequences

Cécile Tran-Kiem<sup>1</sup>, Trevor Bedford <sup>1, 2</sup>

<sup>1</sup> Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA <br>
<sup>2</sup> Howard Hughes Medical Institute, Seattle, WA, USA

## Abstract

Quantifying transmission intensity and heterogeneity is crucial to ascertain the threat posed by infectious diseases and inform the design of interventions. Methods that jointly estimate the reproduction number R and the dispersion parameter k have however mainly remained limited to the analysis of epidemiological clusters or contact tracing data, whose collection often proves difficult. Here, we show that clusters of identical sequences are imprinted by the pathogen offspring distribution, and we derive an analytical formula for the distribution of the size of these clusters. We develop and evaluate a novel inference framework to jointly estimate the reproduction number and the dispersion parameter from the size distribution of clusters of identical sequences. We then illustrate its application across a range of epidemiological situations. Finally, we develop a hypothesis testing framework relying on clusters of identical sequences to determine whether a given pathogen genetic subpopulation is associated with an increased or reduced transmissibility. Our work provides new tools to estimate the reproduction number and transmission heterogeneity from pathogen sequences without requiring building a phylogenetic tree, making it thus easily scalable to large pathogen genome datasets. 

## Overview

This repository contains code and data associated with the above preprint.

### Exploring the size distribution of clusters of identical sequences

In the script **explore_distribution_cluster_identical_sequences.R**, we investigate how epidemiological and evolutionary parameters influence the size distribution of clusters of identical pathogen sequences.

In the script **contribution_largest_transmission_event.R**, we perform a simulation study to understand how much the largest transmission event (within a cluster of identical sequences) contributes to the size of this cluster of identical sequence. 

### Simulation study evaluating our inference framework

**Unconditional on cluster extinction**

**Conditional on cluster extinction**


### Empirical approximation of the probability that transmission occurs before mutation

The script **draw_proba_transmission_before_mutation.R** reproduces the simulation approach developped in the paper to estimate the probability that a transmission event occurs before a mutation event for different pathogens. 

### Application to real data

The following scripts can be used to reproduce the analysis of:
- MERS sequences (**mers.R**)
- Measles sequences from the 2017-2018 Italian outbreak in the Veneto province (**measles.R**)
- SARS-CoV-2 sequences in New-Zealand  (**ncov_NZ.R**)

These scripts only use as input the distribution of clusters of identical sequences. To access the raw data, please refer yourself to the references indicated in the manuscript. To facilitate the application of this method to other datasets, we provide the code developped to generate the size distribution of clusters of identical sequences from an alignment (**get_cluster_distrib_from_fasta.R**). The use of this code is illustrated for an arbitrary fasta file.  

### Simulation study evaluating our transmission advantage statistical framework

### Application to SARS-CoV-2 data in Washington state
