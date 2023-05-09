library(reticulate)
use_python("/home/dnaro/.conda/envs/torchcuda/bin/python")
library(SigProfilerMatrixGeneratorR)

inputfolder = "/g/strcombio/fsupek_data/MMR_BER_Project/Processed_data/Calling/Strelka2/VCFs_SUPEK_24/Discrete_variants"
matrices <- SigProfilerMatrixGeneratorR("Discrete_variants", "GRCh38", inputfolder, plot=T, exome=F,
                                        bed_file=NULL, chrom_based=F, tsb_stat=F, seqInfo=F, cushion=100)
