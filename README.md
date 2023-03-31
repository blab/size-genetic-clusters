# Estimating the reproduction number and transmission heterogeneity from the size distribution of clusters of identical pathogen sequences

Cécile Tran-Kiem<sup>1</sup>, Trevor Bedford <sup>1,2</sup>

<sup>1</sup> Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Center, Seattle, WA, USA <br>
<sup>2</sup> Howard Hughes Medical Institute, Seattle, WA, USA

## Abstract

Quantifying transmission intensity and heterogeneity is crucial to ascertain the threat posed by infectious diseases and inform the design of interventions. Methods that jointly estimate the reproduction number R and the dispersion parameter k have however mainly remained limited to the analysis of epidemiological clusters or contact tracing data, whose collection often proves difficult. Here, we show that clusters of identical sequences are imprinted by the pathogen offspring distribution, and we derive an analytical formula for the distribution of the size of these clusters. We develop and evaluate a novel inference framework to jointly estimate the reproduction number and the dispersion parameter from the size distribution of clusters of identical sequences. We then illustrate its application across a range of epidemiological situations. Finally, we develop a hypothesis testing framework relying on clusters of identical sequences to determine whether a given pathogen genetic subpopulation is associated with an increased or reduced transmissibility. Our work provides new tools to estimate the reproduction number and transmission heterogeneity from pathogen sequences without requiring building a phylogenetic tree, making it thus easily scalable to large pathogen genome datasets. 

## Estimating R and k from user data

To facilitate the application of this method to other datasets, we provide the code developped (i) to generate the size distribution of clusters of identical sequences from an alignment and (ii) to estimate R and k from the generated cluster size distribtion. We illustrate how this may be done starting from an arbitrary FASTA alignment. 

```bash
cd scripts/

## Generate the size distribution of clusters of identical sequences from an alignment. 
# First argument: file path to the user-defined alignment
# Second argument: file path to save the cluster size distribution

Rscript ./get_cluster_distrib_from_fasta.R \
    "../data/synthetic-fasta.fasta" "../results/out_path_clust_size_distrib.rds"

# Note: This script also contains plotting functions to visualize the cluster size distribution.
```
```bash
## Run the inference framework on the generated size distribution of clusters of identical sequences
# First argument: file path of the cluster size distribution
# Second argument: value for the probability p that transmission occurs before mutation
# Third argument: value for the fraction of infections sequenced
# Fourth argument: 1 if confidence intervals are computed. 0 otherwise. This is generally the most computationally intensive part.
# Fifth argument: file path (csv) to save the results of the inference

Rscript ./run_inference_from_distrib.R \
    "../results/out_path_clust_size_distrib.rds" 0.8 0.1 0 "../results/out_path_res_inference.csv"
```

## Overview

This repository contains code and data associated with the above preprint.

### Exploring the size distribution of clusters of identical sequences

In the script [**explore_distribution_cluster_identical_sequences.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/explore_distribution_cluster_identical_sequences.R), we investigate how epidemiological and evolutionary parameters influence the size distribution of clusters of identical pathogen sequences.

In the script [**contribution_largest_transmission_event.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/contribution_largest_transmission_event.R), we perform a simulation study to understand how much the largest transmission event (within a cluster of identical sequences) contributes to the size of this cluster of identical sequence. 

### Simulation study evaluating our inference framework

**Unconditional on cluster extinction**

The script [**sim_study.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/sim_study.R) enables for user-defined values of the reproduction number R, the dispersion parameter k, the proportion of infections sequenced and the probability that transmission occurs before mutation (i) to generate synthetic cluster data and (ii) run our inference framework on these synthetic data. The script [**reproduce_sim_study_figure.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/reproduce_sim_study_figure.R) takes as input a dataframe with the results of the extensive simulation study we performed and produces figures to compare the true and estimated parameter values. 

**Conditional on cluster extinction**
The script [**sim_study_conditional_extinction.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/sim_study_conditional_extinction.R) enables for user-defined values of the reproduction number R, the dispersion parameter k, the proportion of infections sequenced and the probability that transmission occurs before mutation (i) to generate synthetic cluster data and (ii) run our inference framework on these synthetic data conditional on clusters having gone extinct. The script [**reproduce_sim_study_cond_extinction_figure.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/reproduce_sim_study_cond_extinction_figure.R) takes as input a dataframe with the results of the extensive simulation study we performed and produces figures to compare the true and estimated parameter values. 

### Simulation study evaluating our transmission advantage statistical framework
The script [**sim_study_trans_adv.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/sim_study_trans_adv.R) enables for user-defined values of the reproduction number of the non-variant, the transmission advantage and the dispersion parameter (i) to generate synthetic cluster data for both the variant and the non-variant and (ii) run our transmission advantage inference framework on this simulated dataset.

### Empirical estimation of the probability that transmission occurs before mutation

The script [**draw_proba_transmission_before_mutation.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/draw_proba_transmission_before_mutation.R) reproduces the simulation approach developped in the paper to estimate the probability that a transmission event occurs before a mutation event for different pathogens. 

### Application to real data

The following scripts can be used to reproduce the analysis of:
- MERS sequences ([**mers.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/mers.R))
- Measles sequences from the 2017-2018 Italian outbreak in the Veneto province ([**measles.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/measles.R))
- SARS-CoV-2 sequences in New-Zealand  ([**ncov_NZ.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/ncov_NZ.R))
- SARS-CoV-2 variants in Washington state, US ([**ncov_WA.R**](https://github.com/blab/size-genetic-clusters/blob/main/scripts/ncov_WA.R))

These scripts only use as input the distribution of clusters of identical sequences. To access the raw data, please refer yourselves to the references indicated in the manuscript.