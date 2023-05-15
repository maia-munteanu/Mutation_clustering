library(reticulate)
use_python("/usr/bin/python3")
library(SigProfilerExtractorR)

sigprofilerextractor(input_type = "vcf",
                     output = "/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Unclustered_mutations_output_2",
                     input_data = "/g/strcombio/fsupek_cancer1/SV_clusters_project/SigProfiler/Unclustered_mutations", reference_genome="GRCh37", opportunity_genome="GRCh37",
                     minimum_signatures=1, maximum_signatures=20, nmf_replicates=100, cpu=4, gpu=T, batch_size=1,collapse_to_SBS96=F)
