# Notes on Zheng et al PBMC data:
#
#  + Downloaded "Fresh 68k PBMCs (Donor A)" data set from:
#    https://support.10xgenomics.com/single-cell-gene-expression/datasets
#
#  + Associated publication: Zheng et al (2017). Massively parallel
#    digital transcriptional profiling of single cells. Nature
#    Communications 8, 14049. doi:10.1038/ncomms14049
#
#  + Specifically, I downloaded the "Gene/cell matrix (filtered)"
#    tar.gz file. Then I moved the files to
#    data/fresh_68k_pbmc_donor_a_filtered and compressed them with gzip.
#
#  + A summary of the data is given here:
#    http://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/fresh_68k_pbmc_donor_a_web_summary.html
#   
#  + These are the numbers given from that webpage:
#    n = 68,579
#    mean reads per cell = 20,491
#    median read per cell = 1,292
#    median genes per cell = 525
#    number of reads = 1,405,303,609
#
library(Matrix)
library(readr)
