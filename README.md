# flexiMAP
### A regression-based method for discovering differential alternative polyadenylation events in standard RNA-seq data
<br>

## Table of contents

* [Overview](https://github.com/kszkop/flexiMAP#overview)
* [Citing flexiMAP](https://github.com/kszkop/flexiMAP#citing-flexiMAP)
* [License](https://github.com/kszkop/flexiMAP#license)
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

<p align="center">
<img src="https://github.com/kszkop/flexiMAP/blob/master/Figure_s5.png" width="750" />
</p>

------------------------------------------------------------------------

## Citing flexiMAP

Please cite the following article when using __flexiMAP__:

Krzysztof J. Szkop, David S. Moss and Irene Nobeli
***flexiMAP: A regression-based method for discovering differential alternative polyadenylation events in standard RNA-seq data***

<a href="https://doi.org/10.5281/zenodo.3238619"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3238619.svg" alt="DOI"></a>

------------------------------------------------------------------------

## License

[![Packagist](https://img.shields.io/packagist/l/doctrine/orm.svg?maxAge=2592000?style=flat)](https://opensource.org/licenses/MIT)

------------------------------------------------------------------------

## Before starting

### Dependencies

 __flexiMAP__ requires R version >= 3.5.3 and the following packages:

* CRAN
    + dplyr (>= 0.8.1)
    + data.table(>=1.12.2)
    + betareg (>= 3.1-2)
* Bioconductor
    + GenomicRanges (>= 1.36.0)
    + IRanges (>=2.16.0)
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
	
#### InputTable
Dataframe containing names of bam files in the first column, condition in format 1,2,3..(and so on) in the second column. Additional columns with optional covariates for modelling with names of covariates as column names.

###### Format of InputTable dataframe
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
Dataframe containing polyadenylation sites in format: GeneID, chromosome, strand, genomic coordinate. The sites can be predicted by one of the software used for this purpose (reviewed [Szkop et al (2017)](https://onlinelibrary.wiley.com/doi/abs/10.1002/bies.201700090) or can be downloaded from one of the polyadenylation databases (for example: [polyAsite](http://polyasite.unibas.ch), [polyA_DB](http://exon.umdnj.edu/polya_db/), [APASdb](http://genome.bucm.edu.cn/utr/)). Unfortunately, the databases often contain information based on older genome version. Then, we can potentially use one of the software for conversion like CrossMAP(ensembl) available [here](http://crossmap.sourceforge.net) or liftOver(ucsc) available: website [here](https://genome.ucsc.edu/cgi-bin/hgLiftOver) or in Bioconductor, we can use UCSC’s Chain file to apply the liftOver() method provided by package rtracklayer.

###### Format of polyAsite dataframe
	head(polyAsite)
	GeneID Chromosome Strand GenomicCoordinate
	CYP2R1 chr11      +      14910363
	LSM1   chr8       +      38021091
	SIK3   chr11      -      116820211

#### reference
Dataframe containing reference in format: GeneID, transID, chromosome, strand, transcript start, transcript end, CDS start, CDS end

###### Format of reference dataframe
	head(reference)
	GeneID  TransID   Chromosome  Strand  TransStart  TransEnd  CDSstart  CDSend
	BHLHE23 NM_080606 chr20       -       61637330    61638387  61637400  61638126
	CRYBA2  NM_057094 chr2        -       219854911   219858127 219854973 219857898
	FNDC8   NM_017559 chr17       +       33448630    33457751  33448712  33457453


#### TINfilter (optional)
Highly degraded genes with 3'bias shoud be removed from analysis. Thereofre, optionally, a vector of gene symbols to be removed from analysis can be provided to avoid degradation bias. The TIN values (representing degradation rate) can be obtained by running RSeQC software. The genes with low TIN values have high 3' bias and should be removed from analysis. Obviously, it can also be used to remove entries unwanted from other reasons.

	TINfilter = c(CYP2R1, LSM1)

#### num_samples 
minimal number of samples with ratio of (0,1).

<p align="center">
<img src="https://github.com/kszkop/flexiMAP/blob/master/Figure_s1.png" width="750" />
</p>

#### exprFilt20
minimal number counts in short region in all samples, default 20

#### name (optional)
Name of output file if different from default (default: flexiMAP_out)

#### normalise
Default FALSE. Optionally the counts can be normalised for library sizes

### Options connected to the beta-regression
More information about beta-regression can be found [here](https://cran.r-project.org/web/packages/betareg/betareg.pdf)
#### link
default "logit", options: "logit", "probit", "cloglog", "cauchit", "log", "loglog"
#####  type
default "ML", options: "ML", "BC", "BR"
##### link.phi
default "log", options: "identity", "log", "sqrt", see ?betareg for explanation

------------------------------------------------------------------------

## Contacts

krzysztof.szkop@gmail.com

i.nobeli@bbk.ac.uk
