library(reticulate)
use_python("/usr/bin/python3")
library(SigProfilerExtractorR)

sigprofilerextractor(input_type = "vcf",
                     output = "./unclustered_VCFs/",
                     input_data = "./unclustered_VCFs/", reference_genome="GRCh37", opportunity_genome="GRCh37",
                     minimum_signatures=1, maximum_signatures=3, nmf_replicates=100, cpu=1, gpu=F, batch_size=1,make_decomposition_plots=F)
                     
