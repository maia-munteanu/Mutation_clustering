#!/usr/bin/env python3
import os
import sys
from SigProfilerSimulator import SigProfilerSimulator as sigSim
from SigProfilerClusters import SigProfilerClusters as hp

project="Test"
path="./VCFs/"
genome="GRCh37"

def main():
    sigSim.SigProfilerSimulator(project, path, genome, contexts=["96"], simulations=100, chrom_based=True)
    hp.analysis(project, genome, contexts="96", simContext=["96"], input_path=path, subClassify=True, standardVC=True, TCGA=False, sanger=False, probability=True)
if __name__ == '__main__':
        main()


