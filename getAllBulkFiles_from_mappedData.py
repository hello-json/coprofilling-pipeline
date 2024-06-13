#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--mappedSamPath', help='/home/sw2448/palmer_scratch/data/mappedData/GE1224.sam', required=True, type=str)
parser.add_argument('-o', '--outPutDir', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/preProcess/', required=True, type=str)
args = parser.parse_args()
mappedSamPath = args.mappedSamPath
outPutDir = args.outPutDir
os.makedirs(outPutDir, exist_ok=True)
#mappedSamPath = "/home/sw2448/palmer_scratch/data/mappedData/GE1224.sam"
#outPutDir = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/preProcess/"
sampleName = (mappedSamPath.split("/")[-1]).split(".")[0]

cmd1 = "/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools view -bS " + mappedSamPath + " > " + outPutDir + "/" + sampleName + ".bam"
cmd2 = "/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools sort " + outPutDir + "/" + sampleName + ".bam > " + outPutDir + "/sort_" + sampleName + ".bam"
cmd3 = "/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools index " + outPutDir + "/sort_" + sampleName + ".bam"

cmd4 = "/home/sw2448/.conda/envs/theSAMBAMBA/bin/sambamba markdup -r -p -t 2 " + outPutDir + "/sort_" + sampleName + ".bam " + outPutDir + "/sort_" + sampleName + "_clean.bam"

cmd5 = "/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools view -h -o " + outPutDir + "/sort_" + sampleName + "_clean.sam " + outPutDir + "/sort_" + sampleName + "_clean.bam"

os.system(cmd1)
os.system(cmd2)
os.system(cmd3)
os.system(cmd4)
os.system(cmd5)


