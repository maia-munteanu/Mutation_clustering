#!/usr/bin/env python3
import os
import sys
from SigProfilerExtractor import sigpro as sig

def main_function():
      sig.sigProfilerExtractor("vcf", "SBS96", "./VCFs/", reference_genome="GRCh37", opportunity_genome = "GRCh37", context_type = "96", exome = False, 
                               minimum_signatures=1, maximum_signatures=2, nmf_replicates=100, resample = True, batch_size=1, cpu=1)
if __name__=="__main__":
   main_function()
