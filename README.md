# Estimating the reproduction number and transmission heterogeneity from the size distribution of clusters of identical pathogen sequences

Cécile Tran-Kiem<sup>1</sup>, Trevor Bedford <sup>1, 2</sup>

<sup>1</sup> Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA <br>
<sup>2</sup> Howard Hughes Medical Institute, Seattle, WA, USA

## Abstract

xxx

## Overview

This repository contains code and data associated with the above paper.

### Exploring the size distribution of clusters of identical sequences

In the script **explore_distribution_cluster_identical_sequences.R**, we investigate how epidemiological and evolutionary parameters influence the size distribution of clusters of identical pathogen sequences.

In the script **contribution_largest_transmission_event.R**, we perform a simulation study to understand how much the largest transmission event (within a cluster of identical sequences) contributes to the size of this cluster of identical sequence. 

### Simulation study evaluating our statistical framework

### Application to real data

The following scripts can be used to reproduce the analysis of:
- MERS sequences (**mers.R**)
- Measles sequences from the 2017-2018 Italian outbreak in the Veneto province (**measles.R**)
- SARS-CoV-2 sequences in New-Zealand  (**ncov_NZ.R**)

These scripts only use as input the distribution of clusters of identical sequences. To access the raw data, please refer yourself to the references indicated in the manuscript. To facilitate the application of this method to other datasets, we provide the code developped to generate the size distribution of clusters of identical sequences from an alignment (**get_cluster_distrib_from_fasta.R**). The use of this code is illustrated for an arbitrary fasta file.  


### Monitoring changes in variants' transmission characteristics