#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--inDir', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/preProcess/seperatedSam/', required=True, type=str)
#parser.add_argument('-o', '--outDir', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/preProcess/seperatedBam/', required=True, type=str)
args = parser.parse_args()
inDir = args.inDir
#outDir = args.outDir
outDir = inDir.split("seperatedSam")[0] + "/seperatedBam/"
os.system("mkdir " + outDir)

#inDir = "/home/sw2448/palmer_scratch/data/testData/result"
#outDir = "/home/sw2448/palmer_scratch/data/testData/bamResult"

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']

for BarcodeA in oneBarcodeArr:
    for BarcodeB in oneBarcodeArr:
        thisBarcodeName = BarcodeA + "+" + BarcodeB
        thisBarcodeDir = outDir + "/" + thisBarcodeName
        os.system("mkdir " + thisBarcodeDir)
        os.system("/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools view -bS " + inDir + "/" + thisBarcodeName+ ".sam > " + thisBarcodeDir + "/" + thisBarcodeName + ".bam")
        os.system("/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools sort " + thisBarcodeDir + "/" + thisBarcodeName + ".bam > " + thisBarcodeDir + "/sort_" + thisBarcodeName + ".bam")
        os.system("/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools index " + thisBarcodeDir + "/sort_" + thisBarcodeName + ".bam")


