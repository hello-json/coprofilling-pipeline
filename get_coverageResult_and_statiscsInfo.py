#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import subprocess
'''
parser = argparse.ArgumentParser()
parser.add_argument('-w', '--window_len', help='1000000', required=True, type=int)
parser.add_argument('-i', '--input_file_path', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE146/cells/CNVresult1e6/BROWSERFILES/method-edivisive/binsize_1e+06_stepsize_1e+06_CNV.bed.gz', required=True, type=str)
args = parser.parse_args()
window_len = args.window_len
input_file_path = args.input_file_path
'''
input_file_path = "/home/sw2448/palmer_scratch/data/testData/bamResult/"
output_file = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE146/cells/coverageResult/coverage.tsv"
def calculate_average_coverage(bam_file):
    samtools_path = '/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools'
    command_header = [samtools_path, 'view', '-H', bam_file]
    result_header = subprocess.run(command_header, capture_output=True, text=True)
    if result_header.returncode != 0:
        raise Exception("samtools view命令执行失败: " + result_header.stderr)

    chrom_lengths = {}
    for line in result_header.stdout.splitlines():
        if line.startswith('@SQ'):
            parts = line.split('\t')
            chrom = parts[1].split(':')[1]
            length = int(parts[2].split(':')[1])
            chrom_lengths[chrom] = length

    command_depth = [samtools_path, 'depth', bam_file]
    result_depth = subprocess.run(command_depth, capture_output=True, text=True)
    if result_depth.returncode != 0:
        raise Exception("samtools depth命令执行失败: " + result_depth.stderr)

    covered_sites = {}
    for line in result_depth.stdout.strip().split('\n'):
        chrom, _, _ = line.split()[:3]
        if chrom not in covered_sites:
            covered_sites[chrom] = 0
        covered_sites[chrom] += 1

    coverage_info = {}
    total_covered_sites = 0
    total_chrom_length = 0
    for chrom, length in chrom_lengths.items():
        if chrom.startswith('chr'):
            avg_coverage = covered_sites.get(chrom, 0) / length
            coverage_info[chrom] = avg_coverage
            total_covered_sites += covered_sites.get(chrom, 0)
            total_chrom_length += length

    # 只计算以"chr"开头的染色体
    coverage_info['Overall'] = total_covered_sites / total_chrom_length if total_chrom_length else 0

    return coverage_info


#bam_file = "/home/sw2448/palmer_scratch/data/testData/bamResult/ATTGAGGA+CTGAGCCA/sort_ATTGAGGA+CTGAGCCA.bam"
#average_depth_dict = calculate_average_depth(bam_file)
#print(average_depth_dict)

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

oneBarcodeArr = ['AACGTGAT', 'AAACATCG', 'ATGCCTAA', 'AGTGGTCA', 'ACCACTGT', 'ACATTGGC', 'CAGATCTG', 'CATCAAGT', 'CGCTGATC', 'ACAAGCTA', 'CTGTAGCC', 'AGTACAAG', 'AACAACCA', 'AACCGAGA', 'AACGCTTA', 'AAGACGGA', 'AAGGTACA', 'ACACAGAA', 'ACAGCAGA', 'ACCTCCAA', 'ACGCTCGA', 'ACGTATCA', 'ACTATGCA', 'AGAGTCAA', 'AGATCGCA', 'AGCAGGAA', 'AGTCACTA', 'ATCCTGTA', 'ATTGAGGA', 'CAACCACA', 'GACTAGTA', 'CAATGGAA', 'CACTTCGA', 'CAGCGTTA', 'CATACCAA', 'CCAGTTCA', 'CCGAAGTA', 'CCGTGAGA', 'CCTCCTGA', 'CGAACTTA', 'CGACTGGA', 'CGCATACA', 'CTCAATGA', 'CTGAGCCA', 'CTGGCATA', 'GAATCTGA', 'GAAGACTA', 'GAGCTGAA', 'GATAGACA', 'GCCACATA']

all_cells_depthDic_arr = []
all_spots_name = []
for BarcodeA in oneBarcodeArr:
    for BarcodeB in oneBarcodeArr:
        thisBarcodeName = BarcodeA + "+" + BarcodeB
        thisBAM = input_file_path + "/" + thisBarcodeName + "/" + "sort_" + thisBarcodeName + ".bam"
        average_depth_dict = calculate_average_coverage(thisBAM)
        #print(average_depth_dict)
        all_cells_depthDic_arr.append(average_depth_dict)
        all_spots_name.append(str(51 - int(cell_name_dic[BarcodeA])) + "x" + cell_name_dic[BarcodeB])

chrArr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'Overall']
tsvFile = open(output_file, 'ab')
#写表头开始
tsvFile.write("spots".encode())
for one_chrArr_ele in chrArr:
    tsvFile.write(("\t" + one_chrArr_ele).encode())
tsvFile.write(("\n").encode())
#写表头结束
for i in range(len(all_cells_depthDic_arr)):
    tsvFile.write(all_spots_name[i].encode())
    for thisChr in chrArr:
        if thisChr in all_cells_depthDic_arr[i]:
            tsvFile.write(("\t" + str(all_cells_depthDic_arr[i][thisChr])).encode())
        else:
            tsvFile.write(("\t0.0").encode())
    tsvFile.write(("\n").encode())
tsvFile.close()


