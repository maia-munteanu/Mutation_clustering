#!/usr/bin/env python3
import os
import sys
from SigProfilerExtractor import sigpro as sig

#def main_function():
 #     sig.sigProfilerExtractor("vcf", "./Signatures", "./VCFs/", reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "96", exome = False, 
 #                              minimum_signatures=1, maximum_signatures=2, nmf_replicates=100, resample = True, batch_size=1, cpu=1)
#if __name__=="__main__":
#   main_function()

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
def main_function():
     matrices = matGen.SigProfilerMatrixGeneratorFunc("closer", "GRCh37", "./VCFs/",plot=False, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
if __name__=="__main__":
   main_function()
