#!/usr/bin/env python3
import os
import sys
from SigProfilerExtractor import sigpro as sig

project = sys.argv[1]
genome = sys.argv[2]
output = sys.argv[3]

from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
#def main_function():
#     matrices = matGen.SigProfilerMatrixGeneratorFunc("closer", "GRCh37", "./VCFs/",plot=False, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
#if __name__=="__main__":
#   main_function()

def main_function():
     matrices = matGen.SigProfilerMatrixGeneratorFunc(project, genome, output,plot=False, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
if __name__=="__main__":
   main_function()
