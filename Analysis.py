import os
import numpy as np
import matplotlib.pyplot as plt

data_docu = 'real_data_6219'
file_name = 'merged.fa'
data_path = os.path.join(data_docu, file_name)

# value_list = []
# with open(data_path, 'r') as f:
#     lines = f.readlines()
#     for i in range(len(lines)):
#         arr = lines[i].strip()
#         if arr.isnumeric():
#             value_list.append(arr)
#
# file_name = 'merged.fa'
# data_path = os.path.join(data_docu, file_name)
# key_list = []
# with open(data_path, 'r') as f:
#     lines = f.readlines()
#     for line in lines:
#         if (line[0] != '>'):
#             key_list.append(line.strip())
#
# file_name = 'outpout.txt'
# data_path = os.path.join(data_docu, file_name)
# with open(data_path, 'w') as writer:
#     for i in range(len(value_list)):
#         writer.write(key_list[i]+"_"+str(value_list[i])+"\n")

primer_list = set()
with open(data_path, 'r') as f:
    lines = f.readlines()
    for i in range(len(lines)):
        arr = lines[i].strip().split('_')
        if arr[0] == '>oligos':
            primer_pair = tuple([arr[2], arr[3]])
            primer_list.add(primer_pair)

output_docu = 'real_data_6219/Primer'
output_file = 'primer.txt'
output_path = os.path.join(output_docu, output_file)
with open(output_path, 'w') as f:
    for pair in primer_list:
        f.write(pair[0] + ' ' + pair[1] + '\n')





