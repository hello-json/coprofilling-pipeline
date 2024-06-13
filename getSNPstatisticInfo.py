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
outFile = outFileDir + (outFileDir.split("/process_data/")[1]).split("/")[0] + "_snv.tsv"

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']

chromosomesCodingDic = {"chr1": 0, "chr2": 0, "chr3": 0, "chr4": 0, "chr5": 0, "chr6": 0, "chr7": 0, "chr8": 0, "chr9": 0, "chr10": 0, "chr11": 0, "chr12": 0, "chr13": 0, "chr14": 0, "chr15": 0, "chr16": 0, "chr17": 0, "chr18": 0, "chr19": 0, "chr20": 0, "chr21": 0, "chr22": 0, "chrX": 0, "chrY": 0}

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

# 检查input_dir是否包含指定的字符串
contains_special_dir = any(special in resultDir for special in ["GENOHT", "GENOHN", "HGENO", "GE146", "GE1224", "GE838T", "GENO2T", "GENO584", "GENO840", "GE23025", "GE23025new", "GE26355", "GECO146"])

# 定义要输出的染色体列
if contains_special_dir:
    chromosome_columns = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
else:
    chromosome_columns = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX", "chrY"]

all_snp_inAllBarcode = []
tsvFile = open(outFile, 'wb')
tsvFile.write(("cell\t" + "\t".join(chromosome_columns) + "\n").encode())

iii = 0
for BarcodeA in oneBarcodeArr:
    iii += 1
    print(iii)
    for BarcodeB in oneBarcodeArr:
        try:
            chromosomesCodingDic = {chrom: 0 for chrom in chromosome_columns}
            thisBarcodeName = BarcodeA + "+" + BarcodeB
            thisBarcodeCellData = os.path.join(resultDir, thisBarcodeName, thisBarcodeName + ".vcf")
            iniFaObject = open(thisBarcodeCellData)
            try:
                iniFaFile = iniFaObject.read()
            finally:
                iniFaObject.close()
            iniFaFile_rows = iniFaFile.split('\n')
            iniFaFile_rows.pop()
            for one_iniFaFile_row in iniFaFile_rows:
                if not one_iniFaFile_row.startswith("#"):
                    this_chr = one_iniFaFile_row.split("\t")[0]
                    if this_chr in chromosomesCodingDic:
                        #DP = int(((one_iniFaFile_row.split("\t")[7]).split("DP=")[1]).split(";")[0])
                        DP = 1
                        chromosomesCodingDic[this_chr] += DP
            thisCellName = str(51 - int(cell_name_dic[BarcodeA])) + "x" + cell_name_dic[BarcodeB]
            tsvFile.write((thisCellName + "\t").encode())
            for oneChr in chromosome_columns:
                if oneChr == chromosome_columns[-1]:
                    tsvFile.write((str(chromosomesCodingDic[oneChr]) + "\n").encode())
                else:
                    tsvFile.write((str(chromosomesCodingDic[oneChr]) + "\t").encode())
        except:
            thisCellName = str(51 - int(cell_name_dic[BarcodeA])) + "x" + cell_name_dic[BarcodeB]
            print(thisCellName + "      " + BarcodeA + "+" +  BarcodeB)
            tsvFile.write((thisCellName + "\t" + "\t".join(["0"] * len(chromosome_columns)) + "\n").encode())

tsvFile.close()

