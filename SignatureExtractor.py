#!/usr/bin/env python3
import os
import sys
from SigProfilerExtractor import sigpro as sig

def main_function():
      sig.sigProfilerExtractor("matrix", "./Unclustered/Signatures", "./Unclustered/unclustered.SBS96.all", reference_genome="GRCh37", minimum_signatures=1, maximum_signatures=5, nmf_replicates=100, cpu=1)
if __name__=="__main__":
   main_function()

