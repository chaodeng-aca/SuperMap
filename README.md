# SuperMap: Bridging Unpaired Single-Cell Multimodal Data for Integrative Analyses

## Description

**SuperMap** is an R package for the integrative analyses of unpaired single-cell multimodal data. While a substantial amount of single-cell multimodal data has been generated and accumulated, most of these datasets are unpaired, characterized by distinct feature spaces and a lack of cell-wise correspondence. The absence of explicit linkages between modalities poses a fundamental challenge for data integration and interpretation. To address this, we introduce SuperMap, a statistical learning method. SuperMap learns cross-modal mappings to effectively bridge and link different modalities, facilitating a wide range of downstream integrative analysis tasks, including cell-type identification, diagonal integration, regulatory analysis, and trajectory inference.

![SuperMap Figure](SuperMap.png) <!-- Replace with the actual path of the image -->

## Installation

SuperMap is implemented as an R package which can be installed by running the following command in R:

```r
devtools::install_github('chaodeng-aca/SuperMap')
```

## Tutorials

For a step-by-step tutorial, please refer to the `vignette` directory in the repository, which provides usage examples.

* Example 1: [Running on the 10x Multiome Dataset](https://chaodeng-aca.github.io/SuperMap/Running_on_the_10x_Multiome_Dataset.html)
* Example 2: [Running on the ASAP-seq Dataset](https://chaodeng-aca.github.io/SuperMap/Running_on_the_ASAP-seq_Dataset.html)
 
