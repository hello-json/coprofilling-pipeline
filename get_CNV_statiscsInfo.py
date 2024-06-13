#!/home/sw2448/project/conda_envs/thePYTHON/bin/python3.11
# -*- coding: UTF-8 -*-
import os
import gzip
import argparse

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

output_file_path = output_file_dir + "/" + str(window_len) + "_CNV_statiscsInfo.tsv"


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
        ini_cell_name = ((one_input_content_arr_ele.split(" ")[4]).split(".")[0])[-17:]
        #print(ini_cell_name)
        ini_cell_name_arr = ini_cell_name.split("+")
        #print(ini_cell_name_arr)
        cell_name = str(51 - int(cell_name_dic[ini_cell_name_arr[0]])) + "x" + cell_name_dic[ini_cell_name_arr[1]]
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

tsvFile = open(output_file_path, 'ab')
#写表头开始
tsvFile.write("spots".encode())
del input_content_arr[0]
for one_element in input_content_arr:
    if one_element.startswith("track"):
        break
    this_chr = one_element.split("\t")[0]
    for times in range(int((int(one_element.split("\t")[2]) - int(one_element.split("\t")[1])) / window_len)):
        tsvFile.write(("\t" + this_chr).encode())
tsvFile.write(("\n").encode())
#写表头结束
for one_tsv_dic_ele in tsv_dic:
    tsvFile.write(one_tsv_dic_ele.encode())
    for one_one_tsv_dic_ele_value_ele in tsv_dic[one_tsv_dic_ele]:
        tsvFile.write(("\t" + one_one_tsv_dic_ele_value_ele).encode())
    tsvFile.write(("\n").encode())
tsvFile.close()

