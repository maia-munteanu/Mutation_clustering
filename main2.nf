#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2022 IRB Barcelona


params.closer_value = 2000
params.close_value = 10000
params.input_file = "/g/strcombio/fsupek_cancer1/SV_clusters_project/input.csv"
params.output_folder = "/g/strcombio/fsupek_cancer1/SV_clusters_project/New_pipe_results"
params.fasta_ref = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.fasta"
params.hg19 = "/g/strcombio/fsupek_cancer1/SV_clusters_project/hg19.genome"
params.CRG75 = "/home/mmunteanu/reference/CRG75_nochr.bed"

close_bp = params.close_value
closer_bp = params.closer_value
fasta_ref = file(params.fasta_ref)
hg19 = file(params.hg19)
CRG75 = file(params.CRG75)

process serialize_genome {
    publishDir "TMP/"
    conda '/g/strcombio/fsupek_home/dmas/ENV/py36'
    afterScript 'set +u; conda deactivate'
    input:
    file "${params.assembly}.fa" from genome

    output:
    file "${params.assembly}.fa.p" into serial_genome
    file "available_chromosomes.txt" into available_chromosomes

    """
    python -m randommut -M serialize -g ${params.assembly}.fa -a ${params.assembly}

    egrep ">" ${params.assembly}.fa | sed 's/>//' > available_chromosomes.txt
    """
}
