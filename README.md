# flexiMAP
###A regression-based method for discovering differential alternative polyadenylation events in standard RNA-seq data
<br>

## Overview

We present __flexiMAP__ (“flexible Modeling of Alternative PolyAdenylation), a new beta-regression-based method implemented in R, for discovering differential alternative polyadenylation events in standard RNA-seq data. Importantly, __flexiMAP__ allows the modeling of multiple known covariates that often confound the results of RNA-seq data analysis. We show, using simulated data, that __flexiMAP__ is very specific and outperforms in sensitivity existing methods, especially at low fold changes. In addition, the tests on simulated data reveal some hitherto unrecognised caveats of existing methods. 

------------------------------------------------------------------------

## Citing flexiMAP

Please cite the following article when using __flexiMAP__:

Krzysztof J. Szkop, David S. Moss and Irene Nobeli
***flexiMAP: A regression-based method for discovering differential alternative polyadenylation events in standard RNA-seq data***

------------------------------------------------------------------------

## Before starting

### Dependencies

 __flexiMAP__ requires R version >= 3.5.3 and the following packages:

* CRAN
    + dplyr (>= 0.8.1)
    + betareg (>= 3.1-2)
* Bioconductor
    + GenomicRanges (>= 1.36.0)
    + Rsubread (>=1.34.1)
    + Rsamtools (>=2.0.0)

 All dependencies are automatically installed running the code in the next section.

### Installation

 To install __flexiMAP__ directly from GitHub the *devtools* package is required. If not already installed on your system, run
    
    install.packages("devtools")
	
 Otherwise, load _devtools_ and install __flexiMAP__ by
	
	library(devtools)
    install_github("kszkop/flexiMAP", dependencies = TRUE)

------------------------------------------------------------------------

## Contacts

krzysztof.szkop@gmail.com

i.nobeli@bbk.ac.uk
