Installing delimitR
================
Megan L. Smith
2018-05-18

-   [Installing R](#installing-r)
-   [Installing Python](#installing-python)
-   [Installing fastsimcoal2](#installing-fastsimcoal2)
-   [Installing delimitR](#installing-delimitr)
-   [References](#references)

Installing R
============

This package requires R (&gt;=3.5.0). If you do not have R installed, or if you have an earlier version, install R 3.5.0 or later. Visit <https://cran.r-project.org/bin/macosx/> for installation instructions for MacOSx. For installation instructions for Windows, visit <https://cran.r-project.org/bin/windows/base/>.

Installing Python
=================

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
=======================

In order to simulate data, delimitR uses fastsimcoal2 (Excoffier, 2013). To download fastsimcoal2, visit <http://cmpg.unibe.ch/software/fastsimcoal2/>, and follow the installation instructions. delimitR was built using fastsimcoal2 ver 2.6.0.3 - 14.10.17. If you're using a newer version of fastsimcoal2 and encounter issues, please contact me at <megansmth67@gmail.com>.

Installing delimitR
===================

If you have R (&gt;=3.5.0) and python 2.x with sortedcontainers, then you're ready to install delimitR.

To install delimitR, you need to first install Hadley Wikham's devtools package:

    install.packages('devtools', repos = "http://cran.r-project.org")

Download the delimitR source code, and unpack the tar.gz file. Set your working directory to the directory above the delimitR directory that was created when you unpacked the source code.

    devtools::install('delimitR', dependencies = TRUE)

Now, you're ready to run an analysis in delimitR!

References
==========

Excoffier L, Dupanloup I, Huerta-Sanchez E, Sousa VC, Foll M (2013) Robust Demographic Inference from Genomic and SNP Data. PLoS Genet 9(10):e1003905.
