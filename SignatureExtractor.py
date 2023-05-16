#!/usr/bin/env python3
import os
import sys
from SigProfilerExtractor import sigpro as sig

project = sys.argv[1]
matrix = sys.argv[2]
genome = sys.argv[3]
minsig = sys.argv[4]
maxsig = sys.argv[5]

def main_function():
      sig.sigProfilerExtractor("matrix", project, matrix, reference_genome=genome, minimum_signatures=1, maximum_signatures=5, nmf_replicates=100, cpu=8)
if __name__=="__main__":
   main_function()


