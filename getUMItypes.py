#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-

import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--fq2_file', help='input the fq file, eg: /home/sw2448/palmer_scratch/data/testData/GE146_CKDL230016251-1A_HT575DSX5_L1_2.fq', required=True, type=str)
parser.add_argument('-o', '--tsv_file', help='out the tsv file, eg: /home/sw2448/palmer_scratch/data/testData/eachBarcode_UMIamount.tsv', required=True, type=str)
args = parser.parse_args()
fq2path = args.fq2_file
tsvFileName = args.tsv_file

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
barcode_UMI_dic = {}
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
        UMI = line[22:32]
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
            if thisBARCODE not in barcode_UMI_dic:
                barcode_UMI_dic[thisBARCODE] = [UMI]
            else:
                barcode_UMI_dic[thisBARCODE].append(UMI)
    if line_i % 1000000 == 1:
        print(line_i)


cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

tsvFile = open(tsvFileName, 'wb')
tsvFile.write(("cell\tfragments\n").encode())
for one_ele in barcode_UMI_dic:
    thisUMIarr = barcode_UMI_dic[one_ele]
    UMI_amount = len(list(set(thisUMIarr)))
    thisBarcodeName = str(51 - int(cell_name_dic[one_ele.split("+")[0]])) + "x" + cell_name_dic[one_ele.split("+")[1]]
    tsvFile.write((thisBarcodeName + "\t" + str(UMI_amount) + "\n").encode())
tsvFile.close()

