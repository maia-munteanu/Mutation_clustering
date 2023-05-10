#!/usr/bin/env python3
import os
import sys
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

project = sys.argv[1]
genome = sys.argv[2]
output = sys.argv[3]

def main_function():
     matrices = matGen.SigProfilerMatrixGeneratorFunc(project, genome, output, plot=False, exome=False, bed_file=None, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
if __name__=="__main__":
   main_function()
