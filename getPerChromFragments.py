#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import re
import os
import csv
import argparse
from multiprocessing import Pool

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process SAM files and count fragments per chromosome.')
    parser.add_argument('--directory', required=True, help='Directory containing SAM files')
    parser.add_argument('--fq2path', required=True, help='Path to the .fq file')
    return parser.parse_args()

args = parse_arguments()

directory = args.directory
fq2path = args.fq2path

# 如果fq2path中包含“GE146”，则将值设为“GE146”
if "GE146" in fq2path:
    dir_name = "GE146"
else:
    dir_name = directory.split('/')[-3]

output_tsv = dir_name + '.tsv'

allowed_values = ["GENOHT", "GENOHN", "HGENO", "GE146", "GE1224", "GE838T", "GENO2T", "GENO584", "GENO840", "GE23025", "GE23025new", "GE26355", "GECO146"]

if dir_name in allowed_values:
    chromosome_lengths = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
    }
else:
    chromosome_lengths = {
        'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
        'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
        'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
        'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
        'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299, 'chrY': 91744698
    }

oneBarcodeArr = [
    'AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA',
    'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA',
    'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA',
    'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA',
    'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA'
]
readsID_UMI = {}

pattern = re.compile(r'@(.*?)\s')
with open(fq2path, 'r') as fq2:
    for line_i, line in enumerate(fq2, 1):
        if line_i % 4 == 1:
            match = pattern.search(line)
            if match:
                thisKEY = match.group(1)
            else:
                continue
        elif line_i % 4 == 2:
            readsID_UMI[thisKEY] = line[22:32]

def count_fragments_per_chromosome(file_path):
    fragment_counts = {chrom: set() for chrom in chromosome_lengths.keys()}

    with open(file_path, 'r') as file:
        for line in file:
            try:
                if not line.startswith('@'):
                    parts = line.split()
                    read_id, chrom, pos = parts[0], parts[2], int(parts[3])

                    if read_id in readsID_UMI:
                        umi = readsID_UMI[read_id]
                    else:
                        continue

                    fragment_counts[chrom].add(umi)
            except:
                continue

    counts_list = [len(fragment_counts[chrom]) for chrom in chromosome_lengths.keys()]

    return counts_list

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

def process_file(filename):
    if filename.endswith(".sam") and "+" in filename:
        filepath = os.path.join(directory, filename)
        result_array = count_fragments_per_chromosome(filepath)
        file_label = filename[:-4]
        file_label = str(51 - int(cell_name_dic[file_label.split("+")[0]])) + "x" + cell_name_dic[file_label.split("+")[1]]
        return [file_label] + result_array
    return None

filenames = os.listdir(directory)

with Pool() as pool:
    results = pool.map(process_file, filenames)

filtered_results = [result for result in results if result is not None]

# 创建标题行，第一列是spot，后面是染色体序号
header = ['spot'] + list(chromosome_lengths.keys())

# 按照指定顺序排序
sorted_results = []
for i in range(50, 0, -1):
    for j in range(1, 51):
        label = f"{i}x{j}"
        for result in filtered_results:
            if result[0] == label:
                sorted_results.append(result)
                break

# 将标题行和所有结果写入TSV文件
with open(output_tsv, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(header)
    writer.writerows(sorted_results)

