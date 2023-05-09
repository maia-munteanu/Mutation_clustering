library(reticulate)
use_python("/home/dnaro/.conda/envs/torchcuda/bin/python")
library(SigProfilerMatrixGeneratorR)

args=commandArgs(trailingOnly=TRUE)
inputfolder=as.character(args[1])

matrices <- SigProfilerMatrixGeneratorR("Test", "GRCh38", inputfolder, plot=T, exome=F,
                                        bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
