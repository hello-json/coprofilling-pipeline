#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--sampleName', help='GE1224', required=True, type=str)
args = parser.parse_args()
sampleName = args.sampleName

os.system("mkdir /home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/" + sampleName + "/cells")
outdir = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/" + sampleName + "/cells/bamFile/"
os.system("mkdir " + outdir)

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']

for BarcodeA in oneBarcodeArr:
    for BarcodeB in oneBarcodeArr:
        thisBarcodeName = BarcodeA + "+" + BarcodeB
        thisCMD = "cp /home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/" + sampleName + "/preProcess/seperatedBam/" + thisBarcodeName + "/sort_" + thisBarcodeName + ".bam " + outdir
        try:
            os.system(thisCMD)
        except:
            continue

        

