library(reticulate)
use_python("/usr/bin/python3")
library(SigProfilerExtractorR)

sigprofilerextractor(input_type = "matrix",
                     output = "./Unclustered/Signatures",
                     input_data = "./Unclustered", reference_genome="GRCh37", opportunity_genome="GRCh37",
                     minimum_signatures=1, maximum_signatures=3, nmf_replicates=100, cpu=1, gpu=F, batch_size=1,collapse_to_SBS96=F)
