#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona


params.closer_value = 2000
params.close_value = 10000
params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/input.csv"
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/New_pipe_results"
params.reference = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"
params.assembly = "hg19"
params.chr_sizes = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.genome"
params.CRG75 = "/home/mmunteanu/reference/CRG75_nochr.bed"

reference = file(params.reference)
chr_sizes = file(params.chr_sizes)
CRG75 = file(params.CRG75)

process serialize_genome {
    container = '/home/mmunteanu/Randommut.img'
    input: 
    reference
    file "${params.assembly}.fa" from genome

    output:
    file "!{reference}.p" into serial_genome

    """
    randommut -M serialize -g !{reference} -a !{assembly}
    """
}
