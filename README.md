[![Build Status](https://travis-ci.org/rprops/Phenoflow_package.svg?branch=master)](https://travis-ci.org/rprops/Phenoflow_package)
# Phenoflow
R package incorporating many functions for the analysis of microbial flow cytometry data
===============
- **Authors**: Ruben Props [Ruben.Props@UGent.be], Frederiek-Maarten Kerchkof [FrederiekMaarten.Kerckhof@UGent.be]

<p align="justify">The goal of this package is to provide easy functions to analyse flow cytometry data. One of the applications is the derivation of phenotypic diversity estimates, such as depicted below (data available at <a href="http://datadryad.org/resource/doi:10.5061/dryad.m1c04"> doi:10.5061/dryad.m1c04 </a>), which can be compared with 16S rRNA gene amplicon data sets. </p>![alt text][logo] 
[logo]: https://github.com/rprops/PhenoFlow/blob/master/Animation_low_res.gif "Figure 1"


## Available functions

Functions  | Actions
------------| -----------
flowBasis | Function part of the flowFDA package (in development) which performes bivariate kernel density estimations on the phenotypic parameters
Diversity | Calculation of Hill diversities of order 0, 1 and 2 from fingerprint object
Diversity_16S | Calculation of Hill diversities from 16S amplicon data for comparison with <code>Diversity()</code>. 
Evenness | Calculation of pareto evenness (Wittebolle L. et al. (2009)) from fingerprint object (1 = maximum evenness, 0 = minimum evenness)
So | Calculation of Structural Organization parameter (Koch et al. (2014), Frontiers in Microbiology)
CV | Calculation of Coefficient of Variation (CV) of the fingerprint object
beta.diversity.fcm | Non-metric Multidimensional Scaling (NMDS) or PCoA of the phenotypic fingerprints
time.discretization | Function for subsetting .fcs files in time intervals and exporting them as new .fcs files. Designed for the analysis of on-line experiments.
FCS.resample | Resamples sample files from flowSet object to an equal number of cells. Standard is to the minimum sample size.


## Functionalities to be added in the future:
Functions  | Actions
------------| -----------
flowRate | Calculation of growth rates of microbial populations from flow cytometry data.
flowInsilico | 
