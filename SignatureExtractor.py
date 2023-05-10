


sigProfilerExtractor(input_type, out_put, input_data, reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "default", exome = False, 
                         minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, resample = True, batch_size=1, cpu=-1, gpu=False, 
                         nmf_init="random", precision= "single", matrix_normalization= "gmm", seeds= "random", 
                         min_nmf_iterations= 10000, max_nmf_iterations=1000000, nmf_test_conv= 10000, nmf_tolerance= 1e-15, get_all_signature_matrices= False)
