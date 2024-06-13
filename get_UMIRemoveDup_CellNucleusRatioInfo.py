#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import re
import os
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--fq1samPath', help='/home/sw2448/palmer_scratch/data/mappedData/GE13.sam', required=True, type=str)
parser.add_argument('-fq2', '--fq2path', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/GE13_CKDL230024181-1A_HGKY5DSX7_L1_2.fq', required=True, type=str)
#parser.add_argument('-o', '--outDir', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/preProcess/seperatedSam/', required=True, type=str)
args = parser.parse_args()
fq1samPath = args.fq1samPath
fq2path = args.fq2path
#outDir = args.outDir
#outDir = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/" + (fq1samPath.split("/")[-1]).split(".")[0] + "/preProcess/MTUseSeperatedSam/"
#os.system("mkdir " + outDir)
sampleName = (fq1samPath.split("/")[-1]).split(".")[0]
outFileDir = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/" + sampleName + "/cells/cellNucleusRatioResult/"
os.system("mkdir " + outFileDir)
outFile = outFileDir + "cellNucleusRatio.tsv"

'''
def if_UMI_match(a, b):
      #模糊匹配阈值修改！！！！！！
      match_base_amount = 9
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
'''

key_pattern = re.compile(r'@(.*?)\s')
def if_match(a, b, match_base_amount):
    return sum(1 for x, y in zip(a, b) if x == y) >= match_base_amount

#fq2path = "/home/sw2448/palmer_scratch/data/testData/GE146_CKDL230016251-1A_HT575DSX5_L1_2.fq"

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']
readsID_barcode_dic = {}
readsID_UMI_dic = {}

'''
fq2 = open(fq2path, 'rb')
line_i = 0
#barcode_arr = []
UMI_arr = []
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
            UMI = line[22:32]
            flag = 0
            for one_UMI in UMI_arr:
                if if_UMI_match(UMI, one_UMI):
                    flag = 1
                    UMI = one_UMI
                    break
            if flag == 0:
                UMI_arr.append(UMI)
            readsID_UMI_dic[thisKEY] = UMI
            #barcode_arr.append(thisBARCODE)
    if line_i % 1000000 == 1:
        print(line_i)
        sys.stdout.flush()
'''

with open(fq2path, 'r') as fq2:
    line_i = 0
    UMI_arr = set()  # 使用集合来提高查找效率
    while True:
        line = fq2.readline()
        if not line:
            break
        line_i += 1
        if line_i % 4 == 1:
            match = key_pattern.search(line)
            thisKEY = match.group(1)
            #print(thisKEY)
            #sys.stdout.flush()
        if line_i % 4 == 2:
            BARCODE_B = line[32:40]
            BARCODE_A = line[70:78]
            match_flagA = match_flagB = False
            for oneBarcode in oneBarcodeArr:
                if not match_flagB and if_match(oneBarcode, BARCODE_B, 7):
                    match_flagB = True
                    this_BARCODE_B = oneBarcode
                if not match_flagA and if_match(oneBarcode, BARCODE_A, 7):
                    match_flagA = True
                    this_BARCODE_A = oneBarcode
                if match_flagB and match_flagA:
                    break
            if match_flagB and match_flagA:
                thisBARCODE = this_BARCODE_B + "+" + this_BARCODE_A
                readsID_barcode_dic[thisKEY] = thisBARCODE
                UMI = line[22:32]
                if UMI not in UMI_arr:
                    UMI_arr.add(UMI)
                readsID_UMI_dic[thisKEY] = UMI
        #if line_i % 1000000 == 1:
            #print(line_i)
            #sys.stdout.flush()

print("----------------------------------------------------------------------------------------------------------------------------")
sys.stdout.flush()
#print(len(barcode_arr))
#barcode_arr = list(set(barcode_arr))
#print(len(barcode_arr))
#print(readsID_barcode_dic)
#print(len(readsID_barcode_dic))
#sys.stdout.flush()
#获取readsID_barcode字典结束


barcode_chrInfo_dic = {}
barcode_UMIInfo_dic = {}
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
    this_readCHR = line.split("\t")[2]
    try:
        if readsID_barcode_dic[this_readID] not in barcode_chrInfo_dic:
            if (this_readCHR != "chrM") and (this_readCHR != "*") and ("*" not in this_readCHR):
                barcode_chrInfo_dic[readsID_barcode_dic[this_readID]] = {"M": 1, "notM": 0}
            else:
                barcode_chrInfo_dic[readsID_barcode_dic[this_readID]] = {"M": 0, "notM": 1}
            barcode_UMIInfo_dic[readsID_barcode_dic[this_readID]] = [readsID_UMI_dic[this_readID]]
        else:
            if readsID_UMI_dic[this_readID] not in barcode_UMIInfo_dic[readsID_barcode_dic[this_readID]]:
                if (this_readCHR != "chrM") and (this_readCHR != "*") and ("*" not in this_readCHR):
                    barcode_chrInfo_dic[readsID_barcode_dic[this_readID]]["M"] = barcode_chrInfo_dic[readsID_barcode_dic[this_readID]]["M"] + 1
                else:
                    barcode_chrInfo_dic[readsID_barcode_dic[this_readID]]["notM"] = barcode_chrInfo_dic[readsID_barcode_dic[this_readID]]["notM"] + 1
                (barcode_UMIInfo_dic[readsID_barcode_dic[this_readID]]).append(readsID_UMI_dic[this_readID])
    except:
        continue

#print(barcode_chrInfo_dic)
#sys.stdout.flush()

#写入结果文件开始
cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}
tsvFile = open(outFile, 'ab')
tsvFile.write(("cell" + "\t" + "cellNucleusRatio" + "\n").encode())
for oneBarcode1 in oneBarcodeArr:
    for oneBarcode2 in oneBarcodeArr:
        this_cell_key_name = oneBarcode1 + "+" + oneBarcode2
        this_cell_MTratioInfo_dic = barcode_chrInfo_dic[this_cell_key_name]
        this_cell_MTratio = float(this_cell_MTratioInfo_dic['M']) / (float(this_cell_MTratioInfo_dic['M']) + float(this_cell_MTratioInfo_dic['notM']))
        this_cell_barcode_name = str(51 - int(cell_name_dic[oneBarcode1])) + "x" + cell_name_dic[oneBarcode2]
        tsvFile.write((this_cell_barcode_name + "\t" + str(this_cell_MTratio) + "\n").encode())
tsvFile.close()
#写入结果文件结束
