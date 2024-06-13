#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import os
import gzip
import argparse
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-w', '--window_len', help='1000000', required=True, type=int)
parser.add_argument('-i', '--input_file_path', help='/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE146/cells/CNVresult1e6/BROWSERFILES/method-edivisive/binsize_1e+06_stepsize_1e+06_CNV.bed.gz', required=True, type=str)
args = parser.parse_args()
window_len = args.window_len
input_file_path = args.input_file_path

'''
window_len = 100000
window_len = 500000
window_len = 1000000
input_file_path = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE146/cells/CNVresult1e5/BROWSERFILES/method-edivisive/binsize_1e+05_stepsize_1e+05_CNV.bed.gz"
input_file_path = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE146/cells/CNVresult5e5/BROWSERFILES/method-edivisive/binsize_5e+05_stepsize_5e+05_CNV.bed.gz"
input_file_path = "/home/sw2448/palmer_scratch/data/spatialSingleCellEachCloneVariantCluter_project/process_data/GE146/cells/CNVresult1e6/BROWSERFILES/method-edivisive/binsize_1e+06_stepsize_1e+06_CNV.bed.gz"
'''
output_file_dir = input_file_path.split("BROWSERFILES")[0] + "myStatiscsInfoFiles"
os.system("mkdir " + output_file_dir)

output_file_path = output_file_dir + "/" + str(window_len) + "_bulkCNV_statiscsInfo.csv"

sample_name = input_file_path.split("/")[7]

cell_name_dic = {'AACGTGAT': '1', 'AAACATCG': '2', 'ATGCCTAA': '3', 'AGTGGTCA': '4', 'ACCACTGT': '5', 'ACATTGGC': '6', 'CAGATCTG': '7', 'CATCAAGT': '8', 'CGCTGATC': '9', 'ACAAGCTA': '10', 'CTGTAGCC': '11', 'AGTACAAG': '12', 'AACAACCA': '13', 'AACCGAGA': '14', 'AACGCTTA': '15', 'AAGACGGA': '16', 'AAGGTACA': '17', 'ACACAGAA': '18', 'ACAGCAGA': '19', 'ACCTCCAA': '20', 'ACGCTCGA': '21', 'ACGTATCA': '22', 'ACTATGCA': '23', 'AGAGTCAA': '24', 'AGATCGCA': '25', 'AGCAGGAA': '26', 'AGTCACTA': '27', 'ATCCTGTA': '28', 'ATTGAGGA': '29', 'CAACCACA': '30', 'GACTAGTA': '31', 'CAATGGAA': '32', 'CACTTCGA': '33', 'CAGCGTTA': '34', 'CATACCAA': '35', 'CCAGTTCA': '36', 'CCGAAGTA': '37', 'CCGTGAGA': '38', 'CCTCCTGA': '39', 'CGAACTTA': '40', 'CGACTGGA': '41', 'CGCATACA': '42', 'CTCAATGA': '43', 'CTGAGCCA': '44', 'CTGGCATA': '45', 'GAATCTGA': '46', 'GAAGACTA': '47', 'GAGCTGAA': '48', 'GATAGACA': '49', 'GCCACATA': '50'}

tsv_dic = {}
with gzip.open(input_file_path, 'rt') as f:
    input_content = f.read()
#print(input_content)
input_content_arr = input_content.split("\n")
input_content_arr.pop()

for one_input_content_arr_ele in input_content_arr:
    #print(one_input_content_arr_ele)
    if one_input_content_arr_ele.startswith("track"):
        #ini_cell_name = ((one_input_content_arr_ele.split(" ")[4]).split(".")[0])[-17:]
        #print(ini_cell_name)
        #ini_cell_name_arr = ini_cell_name.split("+")
        #print(ini_cell_name_arr)
        #cell_name = str(51 - int(cell_name_dic[ini_cell_name_arr[0]])) + "x" + cell_name_dic[ini_cell_name_arr[1]]
        cell_name = sample_name
        #print(cell_name)
        tsv_dic[cell_name] = []
    else:
        one_input_content_arr_ele_arr = one_input_content_arr_ele.split("\t")
        recycle_num = int((int(one_input_content_arr_ele_arr[2]) - int(one_input_content_arr_ele_arr[1])) / window_len)
        #print(recycle_num)
        copy_num = (one_input_content_arr_ele_arr[3].split("-"))[0]
        for i in range(recycle_num):
            tsv_dic[cell_name].append(copy_num)

#print(tsv_dic)
#for iii in tsv_dic:
#    if len(tsv_dic[iii]) != 6164:
#        print(len(tsv_dic[iii]))

tsvFile = open(output_file_path, 'wb')
#写表头开始
#tsvFile.write("spots".encode())
#del input_content_arr[0]
#for one_element in input_content_arr:
#    if one_element.startswith("track"):
#        break
#    this_chr = one_element.split("\t")[0]
#    for times in range(int((int(one_element.split("\t")[2]) - int(one_element.split("\t")[1])) / window_len)):
#        tsvFile.write(("\t" + this_chr).encode())
#tsvFile.write(("\n").encode())
#tsvFile.write((sample_name + "_" + str(window_len) + "\tcopyNum\n").encode())
tsvFile.write((sample_name + "_" + str(window_len) + ",copyNum\n").encode())
#写表头结束

def generate_array_b(a):
    b = []
    count_dict = {}
    for item in a:
        if item not in count_dict:
            count_dict[item] = 1
            b.append("{}_1".format(item))
        else:
            count_dict[item] += 1
            b.append("{}_{}".format(item, count_dict[item]))
    return b
ini_column1_arr = []
del input_content_arr[0]
for one_element in input_content_arr:
    this_chr = one_element.split("\t")[0]
    for times in range(int((int(one_element.split("\t")[2]) - int(one_element.split("\t")[1])) / window_len)):
        ini_column1_arr.append(this_chr)
#column1_arr = generate_array_b(ini_column1_arr)

column2_arr = tsv_dic[cell_name]
column2_arr = [int(item) for item in column2_arr]


#开始
#threld = 250   #每threld个取一个平均值
threld = 95
def calculate_averages_and_elements(a, b, c):
    # 归类数组a中相同元素对应的数组b中的元素
    grouped_values = defaultdict(list)
    for i, val in enumerate(a):
        grouped_values[val].append(b[i])

    # 计算每组元素的平均值，同时记录每个平均值对应的元素名
    averages = []
    element_names = []

    for key, values in grouped_values.items():
        i = 0
        temp_averages = []  # 临时存储当前元素的平均值
        while i < len(values):
            if i + c <= len(values) or len(values) - i >= c / 2:  # 直接计算平均值
                segment = values[i:i+c]
                avg = sum(segment) / len(segment)
                temp_averages.append(avg)
                i += c
            else:  # 剩余元素少于c的一半
                segment = values[max(0, i - c):i + (len(values) - i)]
                avg = sum(segment) / len(segment)
                if temp_averages:  # 如果之前有计算过平均值，则更新最后一个平均值
                    temp_averages[-1] = avg
                else:
                    temp_averages.append(avg)
                break  # 处理完剩余元素后退出循环

        # 将计算得到的平均值和对应的元素名添加到最终结果中
        averages.extend(temp_averages)
        element_names.extend([key] * len(temp_averages))

    return averages, element_names

averages, element_names = calculate_averages_and_elements(ini_column1_arr, column2_arr, threld)
ini_column1_arr = element_names
column2_arr = averages
column1_arr = generate_array_b(ini_column1_arr)
column2_arr = [str(item) for item in column2_arr]
#结束

#print(column1_arr)
#print(column2_arr)
#print(len(column1_arr))
#print(len(column2_arr))
'''
def min_repetition_count(a):
    # 使用 Counter 统计每个元素出现的次数
    count_dict = Counter(a)
    # 找出最小的重复次数
    min_count = min(count_dict.values())
    return min_count
min_repetition_counts = min_repetition_count(ini_column1_arr)
print(min_repetition_counts)

for one_tsv_dic_ele in tsv_dic:
    tsvFile.write(one_tsv_dic_ele.encode())
    for one_one_tsv_dic_ele_value_ele in tsv_dic[one_tsv_dic_ele]:
        tsvFile.write(("\t" + one_one_tsv_dic_ele_value_ele).encode())
    tsvFile.write(("\n").encode())
tsvFile.close()
'''

for row in range(len(column1_arr)):
    #tsvFile.write((column1_arr[row] + "\t" + column2_arr[row] + "\n").encode())
    tsvFile.write((column1_arr[row] + "," + column2_arr[row] + "\n").encode())
tsvFile.close()

'''
sorted_indices = sorted(range(len(column2_arr)), key=lambda i: column2_arr[i])
column1_arr_sorted = [column1_arr[i] for i in sorted_indices]
column2_arr_sorted = [column2_arr[i] for i in sorted_indices]
plt.scatter(column1_arr_sorted, column2_arr_sorted)
plt.title('Simple Scatter Plot')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.savefig('scatter_plot.pdf')
'''
