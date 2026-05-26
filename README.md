# Umbels mediate faster net diversification rates and upward migration in a Neotropical radiation

This repository contains data and code used to investigate relationships among inflorescence architecture, diversification, and elevation in *Bomarea* (Alstroemeriaceae). Analyses include state-dependent diversification models (HiSSE), correlated trait evolution models, latent liability analyses, and ancestral state estimation analyses. Scripts not used in the final manuscript are labeled with `_not_used`.

## Repository structure

data/
    cleaned trait datasets in .nexus format for each trait
    trimmed phylogenetic trees

data_prep/
    raw trait dataset
    raw phylogenetic tree
    raw *Bomarea* occurrence data downloaded from GBIF
    scripts used to calculate traits and generate .nexus files
    scripts used to integrate trait datasets into a single dataframe
    shared function file used across trait scripts

scripts/

    latent_liability/
        subfolders correspond to the traits included in each analysis
        .xml files used in latent liability analyses
        latent_liability.R and print_latent_liability_functions.R generate .xml files
        latent_liability_figs.R and latent_liability_figs_functions.R generate figures from .log files
        README.md documents manual edits made to final .xml files

    plotting/
        scripts for geographic distribution maps, trait correlation heat maps,
        elevation and inflorescence-type correlations, and HiSSE visualizations

    rb/
        RevBayes scripts for elevation and inflorescence-type correlation analyses,
        HiSSE analyses, and All Rates Different (ARD) ancestral state estimation

## Software
Analyses were done in R 4.4.2, RevBayes, and BEAST.

