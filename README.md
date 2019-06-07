# flexiMAP
### A regression-based method for discovering differential alternative polyadenylation events in standard RNA-seq data
<br>

## Table of contents

* [Overview](https://github.com/kszkop/flexiMAP#overview)
* [Citing flexiMAP](https://github.com/kszkop/flexiMAP#citing-flexiMAP)
* [Before starting](https://github.com/kszkop/flexiMAP#before-starting)
	- [Dependencies](https://github.com/kszkop/flexiMAP#dependencies)
	- [Installation](https://github.com/kszkop/flexiMAP#installation)
	- [Loading](https://github.com/kszkop/flexiMAP#loading)
	- [GettingHelp](https://github.com/kszkop/flexiMAP#getting-help)
* [Usage](https://github.com/kszkop/flexiMAP#usage)
* [Contacts](https://github.com/kszkop/flexiMAP#contacts)

------------------------------------------------------------------------

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

### Loading

 To load __flexiMAP__ run

	library(flexiMAP)

### Getting help

 Next sections illustrate how to make use of __flexiMAP__ by introducing all functions included in the package and reporting most of the data structures and graphical outputs generated with the default options. Similar information are reported in the vignette returned by
 
	browseVignettes("flexiMAP")
 
 For additional examples and further details about the meaning and usage of all parameters in a function run
 
	?function_name

 or

	help(package = flexiMAP)
 
 A complete reference manual is available [here](https://github.com/kszkop/flexiMAP/master/ReferenceManual.pdf).   

 Bugs and errors can be reported at the [issues](https://github.com/kszkop/flexiMAP/issues) page on GitHub. Before filing new issues, please read the documentation and take a look at currently open and already closed discussions.

------------------------------------------------------------------------

## Usage

	flexiMAP(InputTable, path = NULL, covariates = NULL, polyAsite, reference, TINfilter = NULL, num_samples = 1, link = "logit", type = "ML", link.phi = "log", name = "flexiMAP_out", normalise = FALSE, exprFilt = 20)
	
##### InputTable
Dataframe containing names of bam files in the first column, condition in format 1,2,3..(and so on) in the second column. Additional columns with optional covariates for modelling with names of covariates as column names.

##### Format of InputTable dataframe
	head(InputTable)
	Sample         Condition Sex
	sample1.bam    1         1
	sample2.bam    1         2
	sample3.bam    2         1
	sample4.bam    2         2
	
#### path
Path to folder containg indexed bam files

#### covariates
Names of covariates (additional to condition, the same as column names in InputTable)

#### polyAsite
Dataframe containing polyadenylation sites in format: GeneID, chromosome, strand, genomic coordinate. The sites can be predicted by one of the software used for this purpose (reviewed here: https://onlinelibrary.wiley.com/doi/abs/10.1002/bies.201700090) or can be downloaded from one of the polyadenylation databases (for example: [polyAsite](http://polyasite.unibas.ch), [polyA_DB](http://exon.umdnj.edu/polya_db/), [APASdb](http://genome.bucm.edu.cn/utr/)). Unfortunately, the databases often contain information based on older genome version. Then, we can potentially use one of the software for conversion like CrossMAP(ensembl) available [here](http://crossmap.sourceforge.net) or liftOver(ucsc) available: website [here](https://genome.ucsc.edu/cgi-bin/hgLiftOver) or in Bioconductor, we can use UCSC’s Chain file to apply the liftOver() method provided by package rtracklayer.

##### reference

##### TINfilter
##### num_samples 

##### link
#####  type
##### link.phi

##### name
##### normalise
##### exprFilt20



------------------------------------------------------------------------

## Contacts

krzysztof.szkop@gmail.com

i.nobeli@bbk.ac.uk
