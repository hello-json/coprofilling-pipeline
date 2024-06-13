#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11

import pysam
import argparse
import os
import csv
from natsort import natsorted

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_dir', help='input one sample dir containing all their cells SNV result, eg:/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE1224/cells/SNVresult/', required=True, type=str)
parser.add_argument('-t', '--tumor_nontumor_file', help='Path to the tumor_nontumor.tsv file', required=True, type=str)
args = parser.parse_args()
resultDir = args.input_dir
tumor_nontumor_file = args.tumor_nontumor_file

outFileDir = resultDir + "/myStatisticsInfoFiles/"
os.makedirs(outFileDir, exist_ok=True)

# Output file paths
base_filename = (outFileDir.split("/process_data/")[1]).split("/")[0]
outFile_nontumor = os.path.join(outFileDir, f"{base_filename}_non_tumor_PerBinWithChromSnv.tsv")
output_sum_csv_nontumor = os.path.join(outFileDir, f"{base_filename}_non_tumor_PerBinWithChromSnv_sum.csv")
outFile_tumor = os.path.join(outFileDir, f"{base_filename}_tumor_PerBinWithChromSnv.tsv")
output_sum_csv_tumor = os.path.join(outFileDir, f"{base_filename}_tumor_PerBinWithChromSnv_sum.csv")

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']

chromosome_lengths = {
    'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299, 'chrY': 91744698
}

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

# Define the window size
window = 5000000

def count_fragments_per_bin(file_path, window):
    fragment_counts = {chrom: {i: 0 for i in range((length // window) + 1)}
                       for chrom, length in chromosome_lengths.items()}

    with open(file_path, 'r') as file:
        for line in file:
            if not line.startswith('#'):
                parts = line.split()
                chrom, pos, dp_info = parts[0], int(parts[1]), parts[7]
                if chrom in chromosome_lengths:
                    dp = int(dp_info.split('DP=')[1].split(';')[0])
                    bin_index = pos // window
                    fragment_counts[chrom][bin_index] += dp

    counts_list = []
    for chrom in natsorted(chromosome_lengths.keys()):
        counts_list.extend(fragment_counts[chrom].values())

    return counts_list

# Load tumor and non-tumor classification
tumor_nontumor_dict = {}
with open(tumor_nontumor_file, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)  # Skip header
    for row in reader:
        tumor_nontumor_dict[row[0].strip('"')] = row[1]

# Initialize TSV files
filtered_results_nontumor = []
filtered_results_tumor = []
header_written = False
header = ['cell']
for chrom in natsorted(chromosome_lengths.keys()):
    num_bins = (chromosome_lengths[chrom] // window) + 1
    header.extend([f'{chrom}_{i+1}' for i in range(num_bins)])

# Create TSV files for non-tumor and tumor
with open(outFile_nontumor, 'w') as tsvFile_nontumor, open(outFile_tumor, 'w') as tsvFile_tumor:
    for BarcodeA in oneBarcodeArr:
        for BarcodeB in oneBarcodeArr:
            thisBarcodeName = BarcodeA + "+" + BarcodeB
            thisCellName = str(51 - int(cell_name_dic[BarcodeA])) + "x" + cell_name_dic[BarcodeB]
            try:
                thisBarcodeCellData = os.path.join(resultDir, thisBarcodeName, thisBarcodeName + ".vcf")
                result_array = count_fragments_per_bin(thisBarcodeCellData, window)

                if not header_written:
                    tsvFile_nontumor.write('\t'.join(header) + '\n')
                    tsvFile_tumor.write('\t'.join(header) + '\n')
                    header_written = True

                row = [thisCellName] + result_array
                if tumor_nontumor_dict.get(thisCellName, "non-tumor") == "non-tumor":
                    tsvFile_nontumor.write('\t'.join(map(str, row)) + '\n')
                    filtered_results_nontumor.append(row)
                else:
                    tsvFile_tumor.write('\t'.join(map(str, row)) + '\n')
                    filtered_results_tumor.append(row)
            except Exception as e:
                print(f"Error processing {thisBarcodeName}: {e}")
                row = [thisCellName] + ['0'] * (len(header) - 1)
                if tumor_nontumor_dict.get(thisCellName, "non-tumor") == "non-tumor":
                    tsvFile_nontumor.write('\t'.join(row) + '\n')
                    filtered_results_nontumor.append(row)
                else:
                    tsvFile_tumor.write('\t'.join(row) + '\n')
                    filtered_results_tumor.append(row)

# Calculate column sums
def calculate_column_sums(filtered_results, header):
    column_sums = [0] * (len(header) - 1)
    for result in filtered_results:
        for i in range(1, len(result)):
            column_sums[i-1] += int(result[i])
    return column_sums

column_sums_nontumor = calculate_column_sums(filtered_results_nontumor, header)
column_sums_tumor = calculate_column_sums(filtered_results_tumor, header)

# Create sum files for non-tumor and tumor
def create_sum_file(output_sum_csv, header, column_sums):
    with open(output_sum_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['pos', 'snv'])
        for i in range(1, len(header)):
            writer.writerow([header[i], column_sums[i-1]])

create_sum_file(output_sum_csv_nontumor, header, column_sums_nontumor)
create_sum_file(output_sum_csv_tumor, header, column_sums_tumor)

