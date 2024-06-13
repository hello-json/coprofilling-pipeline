#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import re
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--fq1samPath', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/preProcess/sort_GE1224_clean.sam', required=True, type=str)
parser.add_argument('-fq2', '--fq2path', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/GE1224_CKDL230016252-1A_HT575DSX5_L1_2.fq', required=True, type=str)
#parser.add_argument('-o', '--outDir', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/preProcess/seperatedSam/', required=True, type=str)
args = parser.parse_args()
fq1samPath = args.fq1samPath
fq2path = args.fq2path
#outDir = args.outDir
outDir = fq1samPath.split("/sort_")[0] + "/seperatedSam/"
os.system("mkdir " + outDir)


#获取readsID_barcode字典开始

def if_barcode_match(a, b):
      #模糊匹配阈值修改！！！！！！
      match_base_amount = 7
      #a = "ATCGTGAGTC"
      #b = "ATCGTGAAAA"
      flag = 0
      for i in range(len(a)):
          if a[i] == b[i]:
              flag = flag + 1
      if flag >= match_base_amount:
          return(True)
      else:
          return(False)


#fq2path = "/home/sw2448/palmer_scratch/data/testData/GE146_CKDL230016251-1A_HT575DSX5_L1_2.fq"

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']
readsID_barcode_dic = {}
fq2 = open(fq2path, 'rb')
line_i = 0
#barcode_arr = []
while True:
    line = fq2.readline()
    line = line.decode()
    line_i = line_i + 1
    if not line:
        break
    if line_i % 4 == 1:
        #key_pattern = r'@(.*?)\s'
        #match = re.search(key_pattern, line)
        match = re.search(r'@(.*?)\s', line)
        #if match:
        #    thisKEY = match.group(1)
        thisKEY = match.group(1)
    if line_i % 4 == 2:
        BARCODE_B = line[32:40]
        BARCODE_A = line[70:78]
        match_flagA = 0
        match_flagB = 0
        for one_oneBarcodeArr_ele in oneBarcodeArr:
            if if_barcode_match(one_oneBarcodeArr_ele, BARCODE_B):
                match_flagB = 1
                this_BARCODE_B = one_oneBarcodeArr_ele
            if if_barcode_match(one_oneBarcodeArr_ele, BARCODE_A):
                match_flagA = 1
                this_BARCODE_A = one_oneBarcodeArr_ele
            if match_flagB == 1 and match_flagA == 1:
                break
        if match_flagB == 1 and match_flagA == 1:
            thisBARCODE = this_BARCODE_B + "+" + this_BARCODE_A
            readsID_barcode_dic[thisKEY] = thisBARCODE
            #barcode_arr.append(thisBARCODE)
    if line_i % 1000000 == 1:
        print(line_i)
print("----------------------------------------------------------------------------------------------------------------------------")
#print(len(barcode_arr))
#barcode_arr = list(set(barcode_arr))
#print(len(barcode_arr))
#print(len(readsID_barcode_dic))
#获取readsID_barcode字典结束


#切分sam文件开始
print("开始切分sam文件")
#fq1samPath = "/home/sw2448/palmer_scratch/data/testData/sort_GE146_clean.sam"
#outDir = "/home/sw2448/palmer_scratch/data/testData/result"
os.system("/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools view -H " + fq1samPath + " > " + outDir + "/header.sam")

fq1sam = open(fq1samPath, 'rb')
line_i = 0
while True:
    line_i = line_i + 1
    line = fq1sam.readline()
    line = line.decode()
    if not line:
        break
    if line.startswith("@"):
        continue
    this_readID = line.split("\t")[0]
    try:
        thisSamPath = outDir + "/" + readsID_barcode_dic[this_readID] + ".sam"
    except:
        continue
    if os.path.exists(thisSamPath):
        with open(thisSamPath, 'ab') as file:
            file.write(line.encode())
    else:
        os.system("cp " + outDir + "/header.sam " + thisSamPath)
        with open(thisSamPath, 'ab') as file:
            file.write(line.encode())
    if line_i % 100000 == 1:
        print(line_i)
#切分sam文件结束
