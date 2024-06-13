#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import re
import gzip
import os
import csv
from multiprocessing import Pool
from natsort import natsorted

directory = '/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/MBDNA/preProcess/seperatedSam' #文件夹路径最后不用加斜杠
fq2path = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/MBDNA_CKDL240012253-1A_HFF25DSXC_L1_2.fq"
output_tsv = 'MBDNA.tsv'
output_sum_csv = 'MBDNA_sum.csv'
window = 5000000

# readsID_UMI字典开始
oneBarcodeArr = [
    'AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA',
    'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA',
    'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA',
    'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA',
    'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA'
]
readsID_UMI = {}

pattern = re.compile(r'@(.*?)\s')
with open(fq2path, 'r') as fq2:  # 修改文件路径和打开方式
# with gzip.open(fq2path, 'rt') as fq2:
    for line_i, line in enumerate(fq2, 1):
        if line_i % 4 == 1:
            match = pattern.search(line)
            if match:
                thisKEY = match.group(1)
            else:
                continue
        elif line_i % 4 == 2:
            readsID_UMI[thisKEY] = line[22:32]
# readsID_UMI字典结束

#人
'''
chromosome_lengths = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}
'''

#鼠
chromosome_lengths = {
    'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299, 'chrY': 91744698
}


def count_fragments_per_5M_base(file_path):
    # 计算每个染色体的区间数并初始化fragment_counts字典
    fragment_counts = {chrom: {i: set() for i in range((length // window) + 1)}
                       for chrom, length in chromosome_lengths.items()}

    with open(file_path, 'r') as file:
        for line in file:
            try:
                if not line.startswith('@'):
                    parts = line.split()
                    read_id, chrom, pos = parts[0], parts[2], int(parts[3])

                    # 获取UMI值
                    if read_id in readsID_UMI:
                        umi = readsID_UMI[read_id]
                    else:
                        # 如果readsID_UMI字典中不存在该read_id的记录，可以选择跳过或记录错误
                        continue

                    # 计算所在的5M碱基区间索引
                    bin_index = pos // window

                    # 使用UMI来确保独特性
                    fragment_counts[chrom][bin_index].add(umi)
            except:
                continue
    # 转换字典为计数列表
    counts_list = []
    for chrom in natsorted(chromosome_lengths.keys()):
        for bin_index in range(len(fragment_counts[chrom])):
            counts_list.append(len(fragment_counts[chrom][bin_index]))

    return counts_list

#print(count_fragments_per_5M_base("/home/sw2448/palmer_scratch/data/testData/result/CTCAATGA+AAACATCG.sam"))

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

def process_file(filename):
    if filename.endswith(".sam") and "+" in filename:
        filepath = os.path.join(directory, filename)
        result_array = count_fragments_per_5M_base(filepath)
        file_label = filename[:-4]
        file_label = str(51 - int(cell_name_dic[file_label.split("+")[0]])) + "x" + cell_name_dic[file_label.split("+")[1]]
        return [file_label] + result_array
    return None

# 获取目录中的所有文件名
filenames = os.listdir(directory)

# 使用Pool来并行处理文件
with Pool() as pool:
    results = pool.map(process_file, filenames)

# 过滤掉None结果（非有效.sam文件的处理结果）
filtered_results = [result for result in results if result is not None]

# 创建标题行，第一列是spot，后面是带有染色体信息的bin
header = ['spot']
for chrom in natsorted(chromosome_lengths.keys()):
    num_bins = (chromosome_lengths[chrom] // window) + 1
    header.extend([f'{chrom}_{i+1}' for i in range(num_bins)])

# 将标题行和所有结果写入TSV文件
with open(output_tsv, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(header)
    writer.writerows(filtered_results)

# 计算每列的总和
column_sums = [0] * (len(header) - 1)
for result in filtered_results:
    for i in range(1, len(result)):
        column_sums[i-1] += result[i]

# 创建新的sum文件，并将第一行和sum值行转置写入
with open(output_sum_csv, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['pos', 'fragments'])
    for i in range(1, len(header)):
        writer.writerow([header[i], column_sums[i-1]])

