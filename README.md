## R3CPET
3CPET: Finding Co-factor Complexes in Chia-PET experiment using a Hierarchical Dirichlet Process

<a href="http://www.bioconductor.org/packages/devel/bioc/html/R3CPET.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/devel/R3CPET.svg" title="Whether the package is available on all platforms; click for details."></a>
<a href="http://www.bioconductor.org/packages/devel/bioc/html/R3CPET.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/R3CPET.svg" title="How long since the package was first in a released Bioconductor version (or is it in devel only)."></a> <a href="http://bioconductor.org/packages/stats/bioc/R3CPET.html"><img border="0" src="http://www.bioconductor.org/shields/downloads/R3CPET.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/R3CPET/"><img border="0" src="http://www.bioconductor.org/shields/posts/R3CPET.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/devel/bioc/html/R3CPET.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/R3CPET.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

`3CPET` is method based on a non-parametric hierarchical Dirichlet process, to infer the set of protein networks involved in maintaining ChIA-PET chromatin interactions (or any other targeted chromatin interaction data). 3CPET is based on the observation that a protein can be a backbone element in some chromatin loops and dispensable in others. By combining ChIA-PET data with transcription factors binding sites and protein interactions information, 3CPET tries to infer the set of the most probable cofactor complexes involved in maintaining the different spatial interaction contexts.

## Motivations 
The rationale behind 3CPET is that transcription factors can use a distinct combination of coactivators, depending on the genomic and spatial context. We refer to each possible combination of cofactors and their interactions as **Chromatin maintainer Network** (CMN). 

To infer these sets of chromatin maintainer networks (CMNs), first, we build, for each DNA-DNA interaction, a protein-protein interaction (PPI) network connecting the two interacting DNA regions with the help of TF binding site information. Then, we use this set of networks to infer the most enriched coactivator networks. We use a Hierarchical Dirichlet process algorithm to (i) find the number of CMNs and (ii) the proteins involved in each one of them.

## Installation
R3CPET can be downloaded from bioconductor as follow:
```S
## try http if https is not available
source("https://bioconductor.org/biocLite.R")
biocLite("R3CPET")
```

To install the latest development version of R3CPET:
```S
install.packages("devtools")
library("devtools")
install_github("R3CPET","sirusb")
```

## Issues and questions

For issues, you can post them in this git page.
For questions, it is better to post them on [Bioconductor support] (https://support.bioconductor.org/) using the `R3CPET` tag.

### Maintainer: 
Mohamed Nadhir Djekidel :   nde12 at mails dot tsinghua dot edu dot cn
                         or djek.nad at gmail dot com
                         
