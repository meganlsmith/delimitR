delimitR: manual
================
Megan L. Smith
2018-05-18

-   [Overview](#overview)
-   [Installing delimitR](#installing-delimitr)
    -   [Installing R](#installing-r)
    -   [Installing Python](#installing-python)
    -   [Installing fastsimcoal2](#installing-fastsimcoal2)
    -   [Installing delimitR](#installing-delimitr-1)
-   [Input files](#input-files)
    -   [The Observed SFS](#the-observed-sfs)
    -   [The Traits File](#the-traits-file)
-   [How delimitR build models](#how-delimitr-build-models)
    -   [The Guide Tree](#the-guide-tree)
    -   [Other information needed by delimitR](#other-information-needed-by-delimitr)
        -   [Number of species & sample sizes](#number-of-species-sample-sizes)
        -   [Number of SNPs](#number-of-snps)
        -   [Prefix](#prefix)
        -   [Maximum coalescent interval](#maximum-coalescent-interval)
    -   [Setting up the priors](#setting-up-the-priors)
        -   [Population sizes](#population-sizes)
        -   [Divergence times](#divergence-times)
        -   [Migration rates](#migration-rates)
    -   [Setting up the Models](#setting-up-the-models)
-   [Building a Prior](#building-a-prior)
    -   [Simulating data](#simulating-data)
    -   [Binning the mSFS](#binning-the-msfs)
-   [Building a RF Classifier](#building-a-rf-classifier)
-   [Out-of-the-Bag Error Rates](#out-of-the-bag-error-rates)
-   [Selecting a Model](#selecting-a-model)
    -   [Reading and Preparing the Observed Data](#reading-and-preparing-the-observed-data)
    -   [Applying the RF classifier to the observed data.](#applying-the-rf-classifier-to-the-observed-data.)
-   [References](#references)

Overview
========

delimitR is an R-packge for jointly inferring species boundaries and the mode of speciation. delimitR takes as input a multidimensional SFS and compares a model set determined by user specifications to infer under what mode speciation occurred in the focal taxa. delimitR can compare models that include divergence, gene flow, and population size changes.

Installing delimitR
===================

delimitR is written in R, but it relies on fastsimcoal2 and Python to simulate and summarize data.

Installing R
------------

This package requires R (&gt;=3.5.0). If you do not have R installed, or if you have an earlier version, install R 3.5.0 or later. Visit <https://cran.r-project.org/bin/macosx/> for installation instructions for MacOSx. For installation instructions for Windows, visit <https://cran.r-project.org/bin/windows/base/>.

Installing Python
-----------------

Though delimitR is written in R, it depends on Python to generate the binned Site Frequency Spectrum. Python scripts in this package were written in python 2.7.13.

Additionally, the python modules sortedcontainers and collections must be installed. This can be accomplished using pip or conda.

Please note that if you use conda to manage your python it is likely that the python version called from terminal does not match the python version called by a system command in R. To check this, in the terminal run:

    python --version

In R, run:

    system('python --version')

If these match, and if your python is python 2.x, then you should be okay. Otherwise, you'll need to provide a path to the correct python. Make sure this is python 2.x, and make sure that you have the sortedcontainers and collections modules installed. To change the python that R calls, run the code below, replacing "/Users/msmith/anaconda2/bin" with the path to the python version you wish to use.

    Sys.setenv(PATH = paste("/Users/msmith/anaconda2/bin", Sys.getenv("PATH"), sep=":"))

To check that this worked:

    system('python --version')

Installing fastsimcoal2
-----------------------

In order to simulate data, delimitR uses fastsimcoal2 (Excoffier, 2013). To download fastsimcoal2, visit <http://cmpg.unibe.ch/software/fastsimcoal2/>, and follow the installation instructions. delimitR was built using fastsimcoal2 ver 2.6.0.3 - 14.10.17. If you're using a newer version of fastsimcoal2 and encounter issues, please contact me at <megansmth67@gmail.com>.

Installing delimitR
-------------------

If you have R (&gt;=3.5.0) and python 2.x with sortedcontainers, then you're ready to install delimitR.

To install delimitR, you need to first install Hadley Wikham's devtools package:

    install.packages('devtools', repos = "http://cran.r-project.org")

Download the delimitR source code, and unpack the tar.gz file. Set your working directory to the directory above the delimitR directory that was created when you unpacked the source code.

    devtools::install('delimitR', dependencies = TRUE)

Now, you're ready to run an analysis in delimitR!

Input files
===========

delimitR requires two types of input files: a traits file and an observed SFS.

The Observed SFS
----------------

Observed data must be in the format of a multi-dimensional site frequency spectrum (SFS) as implemented in fastsimcoal2 (fsc2; Excoffier, 2013). More information on this format can be found in the fsc2 documentation. An important consideration is population order. Since we need to ensure that our observed SFS is built in the same way as our simulated SFS, we must ensure that population order is consistent between our observed SFS and the priors we set up for the analysis. The following rules must be followed:

1.  Populations must be named with numbers 0 to n-1, where n is the number of populations.

2.  For the purposes of building the observed SFS, populations need to be ordered from population 0 to population n-1.

The file for the observed SFS should have three rows. The first states how many SFS are in the file (in this case one), and gives some information about what the next line contains. The next line has the number of populations followed by the sample size in each population. In this case, we have three populations, and we sampled 10 alleles from each population. Note that this is the number of haploid individuals. If you're working with a diploid species, you will need to treat each haplotype as a separate individual for the purposes of constructing a SFS.

The user should provide the name of the file, omitting the .obs extension. The file should be in the working directory.

For an example of how to construct the observed SFS from pyrad ouptut files, see Example 2 (coming soon).

The observed SFS for this tutorial is named tutorial\_observed\_MSFS.obs, and it is available in the data folder of this repository.

    observedSFS <- 'tutorial_observed_MSFS'

The Traits File
---------------

The next input file is a traits file. This is a tab-delimited file with two-columns. The first line should be a header, and the headings should be 'traits' and 'species'. These species should be named using the integers 0 to n-1, where n is the number of putative species, and arranged in increasing order (0 to n-1). Numbering and ordering is important to ensure that the SFS constructed by fastsimcoal2 for the simulated data matches the one constructed for your observed data. The traits file should be in the working directory, and the user should provide the filename.

The traits file for this tutorial is named tutorial\_traits.txt, and it is available in the data folder of this repository.

    traitsfile <- 'tutorial_traits.txt'

How delimitR build models
=========================

In order to build a set of models, delimitR requires several things from the user, including a guide tree, and priors on divergence times, population sizes, and migration rates.

Once the prior is built under one guide tree and set of priors, the user can change these settings, and create another set of models. As long as models include the same number of putative species and the same assignment of individuals to species, the two model sets can be included in the same analysis downstream. For an example that uses two sets of priors, see Example 3 (coming soon).

Here, we use an example with three populations to set up the analysis.

We need to generate model files for fastsimcoal2. To do this, we use the setup\_fsc2() function. View this documentation and the required input using ?setup\_fsc2. We need to specify several variables here, each of which is described in detail below.

The Guide Tree
--------------

The user must input a guide tree as newick-formatted string. For example, suppose the guide tree is:

<img src="/Users/peglegmeg/Desktop/delimitR/vignettes/myspeciestree.png" width="100 px" />

We specify it as follows:

    observedtree <- '((0,1),2);'

Remember, when specifying guide trees, populations must be named integers from 0 to n-1.

Other information needed by delimitR
------------------------------------

#### Number of species & sample sizes

The user must specify the MAXIMUM number of species that will be tested. This should match the number of tips in the guide tree, and the number of species in the traits file. This number must also be compatible with the observed SFS provided by the user.

    obsspecies<- 3

The user must also provide a vector of sample sizes. Again, this must be compatible with the traits file and the observed SFS. The sample sizes should be specified in order from population 0 to population n-1.

    obssamplesize <- c(10,10,10)

#### Number of SNPs

The user must specify the number of linkage blocks to simulate. For unlinked SNPs, this is equal to the number of SNPs used to build your observed SFS.

    obssnps <- 1500

#### Prefix

The user must also provide a prefix. This will be used to name the fsc2 input files, as well as other output files. This should be unique for all guide tree + prior combinations in a folder, or files will be overwritten.

    obsprefix <- 'tutorial_guidetree1'

#### Maximum coalescent interval

Finally, the user must specify the maximum coalescent interval for which to consider migration. This allows users to narrow their set of models, which becomes essential when more than three potential species are considered. In our three-population example, there are two coalescent intervals. If the user specifies a maximum coalescent interval of 1, then we will test all scenarios that involve secondary contact between population 0 and population 1, leading to four models:

<img src="/Users/peglegmeg/Desktop/delimitR/vignettes/models_max1.png" width="500 px" />

If we specified the maximum coalescent interval to be 2, then we would test models that included gene flow between populations 0 and 1, and between the ancestor of populations 0 and 1 and population 2, leading to a total of seven models:

<img src="/Users/peglegmeg/Desktop/delimitR/vignettes/models_max2.png" width="500 px" />

For the purposes of this tutorial, we will specify a maximum coalescent interval of 1.

    maxinterval <- 1

These are the models that we have defined:

<img src="/Users/peglegmeg/Desktop/delimitR_manual/figures/Models_Ex1.png" width="400 px" />

Setting up the priors
---------------------

#### Population sizes

The user must specify priors on population sizes. These are specified in order from population 0 to population n-1 as a list of vectors. In this case, we specify the priors:

    obspopsizeprior <- list(c(10000,100000),c(10000,100000),c(10000,100000))

The first vector is for population 0, the second for population 1, and the third for population 2. Note that these are in terms of the number of haploid individuals (as specified in the fsc2 documentation).

#### Divergence times

The user must also supply priors for divergence times. Divergence times are given in terms of the number of generations, and must be supplied as a list of vectors.

Order is important when specifying divergence time priors. Consider the tree to be divided into coalescent intervals. For example, the tree below has two intervals (A and B).

<img src="/Users/peglegmeg/Desktop/delimitR/vignettes/coalescentintervals.png" width="300 px" />

Divergence time priors should be provided in order of coalescent interval. First, provide the priors for all divergence events in the first coalescent interval, moving across the tree from left to right. In the figure above, the first prior would be for divergence between species 0 and 1, and this is the only divergence event in the first coalescent interval (interval A in the figure). Then, provide the priors for the second coalescent interval. In our figure, the next divergence time prior would be for divergence between the ancestor of species 0 and 1 and species 2 (coalescent interval B).

By default, delimitR requires that divergence time priors in different coalescent intervals be non-overlapping. This prevents divergence times from ocurring in an order that changes the species tree. If the user desires overlapping confidence intervals, a list of rules specifying the order of divergence times must be provided. For example, if the user wanted overlapping divergence time confidence intervals for the example, they would need to provide the following, where divergence events are labelled Tdiv1$ to Tdivx$ in the order detailed above, where x is the number of divergence events:

    obsdivtimeprior <- list(c(5000,100000),c(50000,1000000))
    myrules <- c('Tdiv2$>Tdiv1$')

For the purpose of the tutorial, we will use non-overlapping intervals.

    obsdivtimeprior <- list(c(50000,100000),c(500000,1000000))

#### Migration rates

The user must also specify a prior on migration rates. Currently, the program only allows one prior for all migration rates in the default model sets.

    obsmigrateprior <- list(c(0.000005,0.00005))

Setting up the Models
---------------------

Now, we are ready to generate the .tpl and .est files that describe these models.

    library(delimitR)
    setup_fsc2(tree=observedtree,
               nspec=obsspecies,
               samplesizes=obssamplesize,
               nsnps=obssnps,
               prefix=obsprefix,
               maxinterval=maxinterval,
               popsizeprior=obspopsizeprior,
               divtimeprior=obsdivtimeprior,
               migrateprior=obsmigrateprior)

Building a Prior
================

Now that we have our models set up, we are ready to simulate data, bin the SFS, and create a prior!

Simulating data
---------------

Now, we are ready to simulate data under each of our models. For this, we use the fastsimcoalsims() function. This function requires as input the prefix used to generate teh model files, the path to fastsimcoal2, and the number of replicates we wish to generate under each model. Generally, a minimum of 10,000 replicates under each model should be simulated.

    fastsimcoalsims(prefix=obsprefix,
                    pathtofsc='../fsc26',
                    nreps=10000)

Binning the mSFS
----------------

When all polymorphic sites are unlinked and biallelic (Gutenkunst, et. al. 2009), the mSFS is a complete summary of the data. However, inferences based on the mSFS may be inaccurate when too few segregating sites are sampled (Terhost & Song, 2015). To address this problem, we used a binning strategy to further summarize the mSFS following Smith et al. (2017) to generate the binned SFS (bSFS).

The user must specify how many bins will be used to summarize the SFS. This number should not be greater than the sample size of the population with the fewest samples, as this results in sparse sampling of the SFS. Large values lead to a more complete summary of the data, but also lead to a more sparsely sampled SFS and increased computation times. Users should use exploratory analyses to find the optimal trade-off between computational time and error rates. The user must specify the number of classes that will be used in the binned SFS.

    nclasses <- 5

Now, we are ready to create the prior using the function makeprior(). The user must provide the prefix used to name the model files, the number of species, the number of classes to be included in the SFS, a path to the working directory, the name of the traits file, the threshold, the name of the folder to store the prior in, and the number of cores to use.

The threshold is used when the observed mSFS is built using a downsampling approach. For example, in the Example 2 (coming soon), we use a threshold of 50%. This means that only SNPs that are sequenced in at least 50% of alleles in each population are used to build the mSFS. SNPs that are sequenced in more than 50% of the individuals are randomly downsampled. The threshold provided here should be the same as the threshold used to perform the downsampling when the mSFS was constructed. If no downsampling was used, the threshold should be set to 100.

The user must provide a folder where the prior will be stored. This folder should not exist prior to running this code, or should be empty.

    FullPrior <- makeprior(prefix=obsprefix,
                           nspec=obsspecies,
                           nclasses=nclasses,
                           getwd(),
                           traitsfile = traitsfile,
                           threshold=100, 
                           thefolder = 'Prior',
                           ncores = 2)

This should create a Binned\*.obs file for each of the models evaluated, and should return a dataframe containing the prior that will be used to build the RF classifier. At this stage, we have a lot of extraneous files and folders in our working directory. By default, these are not removed by delimitR, because some users may want them for downstream analyses. However, these files can take up a lot of space. If you do not need them, the function clean\_working() can be used to remove these from the working directory.

    clean_working(prefix=obsprefix)

There's one final step to cleaning up our Prior. We want to remove rows that have zero variance. For example, if no SNPs are ever observed in a certain bin of the SFS, across all models and all simulations, this bin adds nothing to our analysis, and we must remove it before building the RF classifier. We do this using the function Prior\_reduced().

    ReducedPrior <- Prior_reduced(FullPrior)

Building a RF Classifier
========================

We use the simulated data to construct a Random Forest (RF) classifier, in which the bins of the bSFS are the predictor variables and the model used to simulate the data is the response variable. To build the RF classifier, delimitR uses the R package “abcrf” (Pudlo et al, 2015). The RF classifier consists of a user-defined number N of decision trees. Each decision tree is constructed from a subset of the prior, and at each node in each decision tree a bin of the bSFS is considered, and a binary decision rule is constructed based on the number of SNPs in the bin. When this classifier is applied to new datasets, we move down nodes until we reach the leaves of the tree, which, in this case, are model indices. Each decision tree votes for a model, and the model receiving the largest portion of votes is selected as the best model.

    myRF <- RF_build_abcrf(ReducedPrior,FullPrior,500)
    myRF
    plot(myRF, training = ReducedPrior)

This returns an RF object and prints a summary of this object, that includes overall oob error rates, and error rates for each model. We'll discuss these in detail in the next section.

Out-of-the-Bag Error Rates
==========================

We use out-of-the-bag (oob) error rates to assess the power of the RF classifier. Since only a portion of the prior is used for the construction of each decision tree, we can take an element of the prior, consider only decision trees constructed without reference to that element, and calculate how often we choose an incorrect model.

The out-of-the-bag error rates from this tutorial are reported in the table below. Note that these may differ slightly when you run this example, because different RF classifiers built using the same data will differ slightly. Below, error rates are reported as proportions.

|    Model| error rate |
|--------:|:-----------|
|  overall| 0.0085     |
|        1| 0.0000     |
|        2| 0.0004     |
|        3| 0.0105     |
|        4| 0.0231     |

Overall error rates are very low for this example, as are error rates for each model. The error rate is highest for model 4 (the model that includes secondary contact). We can interpret this error rate to mean that, when data are simulated under model 4, 2.31 % of the time, a model other than model 4 is selected as the best model. We can break this down further and see which models are chosen when errors occur. Here, we report the number of datasets that were classified as belonging to each model.

|  Generating Model| Model 1 | Model 2 | Model 3 | Model 4 |
|-----------------:|:-------:|:-------:|:-------:|:--------|
|           Model 1|  10000  |    0    |    0    | 0       |
|           Model 2|    0    |   9996  |    0    | 4       |
|           Model 3|    0    |    0    |   9895  | 105     |
|           Model 4|    0    |    14   |   217   | 9769    |

For all 10,000 datasets generated under Model 1, Model 1 was selected as the best model. Of the 10,000 datasets generated under Model 4, Model 4 was chosen as the best model for 9,769, Model 3 was chosen for 217, and Model 2 was chosen for 14.

Selecting a Model
=================

Reading and Preparing the Observed Data
---------------------------------------

We need to get the observed data into the correct format. For this, we need the function prepobserved(). The user must provide:

1.  The name of the observed SFS.

2.  The Full Prior (created with the function makeprior()).

3.  The Reduced Prior (created with the function Prior\_reduced()).

4.  The number of classes to use for binning (should match that used to build the Full Prior).

5.  The number of putative species.

6.  The name of the traits file.

7.  The threshold used for downsampling.

<!-- -->

    myobserved <- prepobserved(
      observedSFS,
      FullPrior,
      ReducedPrior,
      nclasses,
      obsspecies,
      traitsfile=traitsfile,
      threshold = 100)

Applying the RF classifier to the observed data.
------------------------------------------------

Now, we're ready to apply the RF classifier to the observed data. To do this, we use the function RF\_predict\_abcrf(), which requires the RF object, the observed dataset, the Reduced Prior, the Full Prior, and the number of trees, which should match that used to construct the classifier.

    prediction <- RF_predict_abcrf(myRF, myobserved, ReducedPrior, FullPrior, 500)
    prediction

This returns the proportion of trees in the RF that voted for each model.

|  Model| prop. votes |
|------:|:------------|
|      1| 0.000       |
|      2| 0.006       |
|      3| 0.002       |
|      4| 0.992       |

We also regress over out of the bag error rates, following (Pudlo et al., 2015), to estimate the posterior probability of the best model, which is &gt; 0.999 in this case, giving us high confidence in our results. For this tutorial, the observed data was simulated under Model 4, so we know that we choose the correct model!

References
==========

Excoffier L, Dupanloup I, Huerta-Sanchez E, Sousa VC, Foll M (2013) Robust Demographic Inference from Genomic and SNP Data. PLoS Genet 9(10):e1003905.

Gutenkunst RN, Hernandez RD, Williamson SH, Bustamante CD (2009) Inferring the Joint Demographic History of Multiple Populations from Multidimensional SNP Frequency Data. PLoS Genet 5(10):e1000695.

Terhorst J, Song YS (2015) Fundamental limits on the accuracy of demographic inference based on the sample frequency spectrum. Proc Natl Acad Sci 112(25):7677–7682.

Smith ML, Ruffley M, Tank DC, Sullivan J, Carstens BC (2017) Demographic model selection using random forests and the site frequency spectrum. 26(17):4562–4573.

Pudlo P, et al. (2015) Reliable ABC model choice via random forests. Bioinformatics 32(6):859–866.
