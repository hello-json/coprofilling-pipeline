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
output_file = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE146/cells/depthResult/depth.tsv"
def calculate_average_depth(bam_file):
    # 使用指定的samtools路径
    samtools_path = '/home/sw2448/project/conda_envs/theSAMTOOLS/bin/samtools'
    command = [samtools_path, 'depth', bam_file]
    result = subprocess.run(command, capture_output=True, text=True)

    # 检查命令是否成功执行
    if result.returncode != 0:
        raise Exception("samtools depth命令执行失败: " + result.stderr)

    # 解析samtools depth的输出
    depth_data = {}
    total_depth = 0
    total_count = 0
    for line in result.stdout.strip().split('\n'):
        chrom, pos, depth = line.split()
        depth = int(depth)
        if chrom not in depth_data:
            depth_data[chrom] = {'total_depth': 0, 'count': 0}
        depth_data[chrom]['total_depth'] += depth
        depth_data[chrom]['count'] += 1
        total_depth += depth
        total_count += 1

    # 计算每条染色体的平均深度，并存储到字典中
    average_depth_dict = {}
    for chrom in depth_data:
        avg_depth = depth_data[chrom]['total_depth'] / depth_data[chrom]['count']
        average_depth_dict[chrom] = avg_depth

    # 计算整个样本的平均测序深度并添加到字典中
    if total_count > 0:
        average_depth_dict['Overall'] = total_depth / total_count
    else:
        average_depth_dict['Overall'] = 0

    return average_depth_dict

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
        average_depth_dict = calculate_average_depth(thisBAM)
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


