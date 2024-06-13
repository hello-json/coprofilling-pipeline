#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11

import pysam
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', help='input one sample dir containing all their cells SNV result, eg:/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/cells/SNVresult/', required=True, type=str)
args = parser.parse_args()
resultDir = args.input_dir

outFileDir = resultDir + "/myStatiscsInfoFiles/"
os.system("mkdir " + outFileDir)
outFile = outFileDir + (outFileDir.split("/process_data/")[1]).split("/")[0] + "_snvEachVariantReducedAccuracy1.tsv"
#resultDir = "/home/sw2448/palmer_scratch/data/testData/freebayesResult"

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']


all_cellInfoDic_arr = []
all_cellInfoDic_keyArr = []
for BarcodeA in oneBarcodeArr:
    for BarcodeB in oneBarcodeArr:
        thisCellInfoDic = {}
        thisBarcodeName = BarcodeA + "+" + BarcodeB
        thisBarcodeCellData = resultDir + "/" + thisBarcodeName + "/" + thisBarcodeName + ".vcf"
        iniFaObject = open(thisBarcodeCellData)
        try:
            iniFaFile = iniFaObject.read( )
        finally:
            iniFaObject.close( )
        iniFaFile_rows = iniFaFile.split('\n')
        iniFaFile_rows.pop()
        #print (iniFaFile_rows)
        for one_iniFaFile_row in iniFaFile_rows:
            if not one_iniFaFile_row.startswith("#"):
                this_chr = one_iniFaFile_row.split("\t")[0]
                if "." not in this_chr:
                    this_pos = one_iniFaFile_row.split("\t")[1]
                    this_ref = one_iniFaFile_row.split("\t")[3]
                    this_alt = one_iniFaFile_row.split("\t")[4]
                    this_key = this_chr + "_" + this_pos
                    #DP = int(((one_iniFaFile_row.split("\t")[7]).split("DP=")[1]).split(";")[0])
                    DP = 1
                    thisCellInfoDic[this_key] = DP
                    all_cellInfoDic_keyArr.append(this_key)
        all_cellInfoDic_arr.append(thisCellInfoDic)

#print(len(all_cellInfoDic_arr))
#print(len(all_cellInfoDic_keyArr))
all_cellInfoDic_keyArr = list(set(all_cellInfoDic_keyArr))
#print(len(all_cellInfoDic_keyArr))

#begin to sort the all_cellInfoDic_keyArr
b = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"]
def custom_sort(element):
    parts = element.split('_')
    return b.index(parts[0]), int(parts[1])
sorted_all_cellInfoDic_keyArr = sorted(all_cellInfoDic_keyArr, key=custom_sort)
#print(sorted_all_cellInfoDic_keyArr)
#end to sort the all_cellInfoDic_keyArr

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

tsvFile = open(outFile, 'ab')
tsvFile.write(("cell").encode())
for one_sorted_all_cellInfoDic_keyArr_ele in sorted_all_cellInfoDic_keyArr:
    tsvFile.write(("\t" + one_sorted_all_cellInfoDic_keyArr_ele).encode())
tsvFile.write(("\n").encode())

i = -1
for BarcodeA in oneBarcodeArr:
    for BarcodeB in oneBarcodeArr:
        i = i + 1
        print(i)
        #thisBarcodeName = BarcodeA + "+" + BarcodeB
        thisBarcodeName = str(51 - int(cell_name_dic[BarcodeA])) + "x" + cell_name_dic[BarcodeB]
        tsvFile.write((thisBarcodeName).encode())
        for one_sorted_all_cellInfoDic_keyArr_ele in sorted_all_cellInfoDic_keyArr:
            if one_sorted_all_cellInfoDic_keyArr_ele in all_cellInfoDic_arr[i]:
                this_value = str(all_cellInfoDic_arr[i][one_sorted_all_cellInfoDic_keyArr_ele])
            else:
                this_value = "0"
            tsvFile.write(("\t" + this_value).encode())
        tsvFile.write(("\n").encode())
tsvFile.close()

