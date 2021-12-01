import os
import numpy as np
import matplotlib.pyplot as plt

data_docu = './'
file_name = 'combined_clusters_with_references.txt'
data_path = os.path.join(data_docu, file_name)

value_list = []
with open(data_path, 'r') as f:
    lines = f.readlines()
    for i in range(len(lines)):
        arr = lines[i].strip()
        if arr.isnumeric():
            value_list.append(arr)

file_name = 'merged.fa'
data_path = os.path.join(data_docu, file_name)
key_list = []
with open(data_path, 'r') as f:
    lines = f.readlines()
    for line in lines:
        if (line[0] != '>'):
            key_list.append(line.strip())
    
file_name = 'outpout.txt'
data_path = os.path.join(data_docu, file_name)
with open(data_path, 'w') as writer:
    for i in range(len(value_list)):
        writer.write(key_list[i]+"_"+str(value_list[i])+"\n")







    
