# Molecular Evolution Analyses

In the context of my Ph.D., I am studying the molecular evolution of the mud snail genus *Ecrobia*. Five *Ecrobia* species are included in this study: 
1. *Ecrobia grimmi* (Ecgr)
2. *Ecrobia maritima* (Ecma)
3. *Ecrobia spalatiana* (Ecspa)
4. *Ecrobia ventrosa* (Ecve)
5. *Ecrobia truncata* (Ectr)

And two outgroups:
1. *Peringia ulvae* (Peul)
2. *Salenthydrobia ferrerii* (Safe)


## Data

My data include two mitochondrial markers:
1. Cytochrome c oxidase subunit I (COI): 359 individuals
2. 16S ribosomal RNA (r16S): 74 individuals

And one nuclear marker:

3. Internal transcribed spacer 2 (ITS2): 46 individuals

## Prerequisites

To estimate maximum likelihood trees with RAxML, this program needs to be previously installed, and the path of the executable needs to be given.

## Installing

I set the working directory to source file location.

## Steps

The whole script is still under progress. The following steps are already included:
1. Import DNA sequences (including outgroups)
2. Exploratory data analyses
  * Base frequency
  * Base frequency for each individual
  * Base frequency between genes for all species pooled - WIP
3. Multiple sequence alignment
4. Best-fit substitution models estimation
5. Saturation test
6. Merge the alignments
7. Infer trees
  * Brute force
  * Evolutionary distances
  * Distance-based phylogenetic methods
  * Maximum parsimony
  * Maximum likelihood methods
  * Bayesian inference
  * Codon model - WIP
  * The class 'phylo' - WIP
  * Visualization methods - WIP
  * Plot trees - WIP

The following steps are not included yet:

8. Tree comparison
9. Bootstrap
10. Divergence time estimation
11. Coalescence analysis
12. List trees
13. Tree space
14. phytools

## Written and Compiled With

* Sublime Text 2 version 2.0.2
* R version 3.4.1
* RStudio version 1.1.419

## Author

* **Justine Vandendorpe**

## Funding

My project is funded by the European Union's Horizon 2020 research and innovation programme under the Marie Sklodowska-Curie grant agreement No 642973.

## Acknowledgments

This script was started during the course Phylogenetic Analysis Using R organized by Transmitting Science in March 2017. I would like to thank Emmanuel Paradis and Klaus Schliep for teaching during this course. I also found a lot of information in Emmanuel Paradis's book 'Analysis of Phylogenetics and Evolution with R' (Second Edition). 