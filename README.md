[![Build Status](https://travis-ci.org/rprops/Phenoflow_package.svg?branch=master)](https://travis-ci.org/rprops/Phenoflow_package)
# Phenoflow
*******************
- **Author**: Ruben Props [Ruben.Props@UGent.be]
- **Contributor**: Frederiek-Maarten Kerckhof [FrederiekMaarten.Kerckhof@UGent.be]

The goal of this package is to provide functions for advanced analysis of microbial flow cytometry data. One of the applications is the derivation of phenotypic diversity estimates, such as depicted below (data available at <a href="http://datadryad.org/resource/doi:10.5061/dryad.m1c04"> doi:10.5061/dryad.m1c04</a>), which can be compared with 16S rRNA gene amplicon data sets. 

*******************
- **Tutorial can be found [here](https://github.com/rprops/Phenoflow_package/wiki/1.-Phenotypic-diversity-analysis).**
*******************

If you use this package, please consider citing the original publication:  

Props R, Monsieurs P, Mysara M, Clement L, Boon N (2016). **Measuring the biodiversity of microbial communities by flow cytometry**. Methods in Ecology and Evolution 7: 1376-1385.

Install the package:
```R
library("devtools")
install_github("CMET-UGent/Phenoflow_package")
```
*******************

![alt text][logo]

[logo]: https://github.com/rprops/PhenoFlow/blob/master/Animation_low_res.gif "Figure 1"

## Available functions

Functions  | Actions
------------| -----------
flowBasis | Function part of the flowFDA package (in development) which performes bivariate kernel density estimations on the phenotypic parameters
Diversity | Calculation of Hill diversities of order 0, 1 and 2 from fingerprint object
Diversity_rf | More accurate (i.e., smaller errors) calculation of Hill diversities from flowSet. It is slower than Diversity().
Diversity_16S | Calculation of Hill diversities from 16S amplicon data for comparison with <code>Diversity()</code>. 
Evenness | Calculation of pareto evenness (Wittebolle L. et al. (2009)) from fingerprint object (1 = maximum evenness, 0 = minimum evenness)
So | Calculation of Structural Organization parameter (Koch et al. (2014), Frontiers in Microbiology)
CV | Calculation of Coefficient of Variation (CV) of the fingerprint object
beta_diversity_fcm | Non-metric Multidimensional Scaling (NMDS) or PCoA of the phenotypic fingerprints
time_discretization | Function for subsetting .fcs files in time intervals and exporting them as new .fcs files. Designed for the analysis of on-line experiments.
FCS_resample | Resamples sample files from flowSet object to an equal number of cells. Standard is to the minimum sample size.
FCS_clean | Denoises all samples by means of [automated identification and removal of fluorescence anomalies in flow cytometry data](http://onlinelibrary.wiley.com/doi/10.1002/cyto.a.22837/full)

## Functionalities to be added in the future:
Functions  | Actions
------------| -----------
flowRate | Calculation of growth rates of microbial populations from flow cytometry data.
