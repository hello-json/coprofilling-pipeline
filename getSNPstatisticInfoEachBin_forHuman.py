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
outFile = outFileDir + (outFileDir.split("/process_data/")[1]).split("/")[0] + "_PerBinSnv.tsv"

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']

chromosome_lengths = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

def count_fragments_per_bin(file_path, window=50000000):
    # 计算每个染色体的区间数并初始化fragment_counts字典
    fragment_counts = {chrom: {i: 0 for i in range((length // window) + 1)}
                       for chrom, length in chromosome_lengths.items()}

    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.split()
                chrom, pos, dp_info = parts[0], int(parts[1]), parts[7]
                if chrom in chromosome_lengths:
                    dp = int(dp_info.split('DP=')[1].split(';')[0])
                    #dp = 1
                    bin_index = pos // window
                    fragment_counts[chrom][bin_index] += dp
    
    # 转换字典为计数列表
    counts_list = []
    for chrom in sorted(chromosome_lengths.keys()):
        counts_list.extend(fragment_counts[chrom].values())
    
    return counts_list

# 获取所有文件名
filenames = os.listdir(resultDir)

# 初始化tsv文件
with open(outFile, 'w') as tsvFile:
    header_written = False
    
    for BarcodeA in oneBarcodeArr:
        for BarcodeB in oneBarcodeArr:
            try:
                thisBarcodeName = BarcodeA + "+" + BarcodeB
                thisBarcodeCellData = os.path.join(resultDir, thisBarcodeName, thisBarcodeName + ".vcf")

                # 统计每个bin的fragments
                result_array = count_fragments_per_bin(thisBarcodeCellData)
                
                # 写header行
                if not header_written:
                    header = ['cell'] + [f'bin{i+1}' for i in range(len(result_array))]
                    tsvFile.write('\t'.join(header) + '\n')
                    header_written = True
                
                # 写结果行
                thisCellName = str(51 - int(cell_name_dic[BarcodeA])) + "x" + cell_name_dic[BarcodeB]
                tsvFile.write(thisCellName + '\t' + '\t'.join(map(str, result_array)) + '\n')
            except Exception as e:
                print(f"Error processing {thisBarcodeName}: {e}")
                tsvFile.write(f"{thisCellName}\t" + '\t'.join(['0'] * len(result_array)) + '\n')

