from Utils import *
import os
import numpy as np
import matplotlib.pyplot as plt

docu = 'real_data_6219'
primer_path = os.path.join(docu, 'Primer', 'primer.txt')
primers = read_primers(primer_path)
reference = os.path.join(docu, 'raw_data', 'outpout.txt')
with open(reference, 'r') as f: ref_str = f.readlines()

result_path = os.path.join(docu, 'Output', 'output-id.csv')
result = []
distance = []
wrong_distance = []
with open(result_path, 'r') as f:
    lines = f.readlines()
    for line in lines:
        arr = line.strip().split(',')
        result.append(int(arr[0]))
        distance.append(int(arr[1]))

index_path = os.path.join(docu, 'dna_pool-index.csv')
with open(index_path, 'r') as f: index = f.readlines()
index = [int(a.strip()) for a in index]

acc = []
miss = 0
for i in range(len(result)):
    id = index[i]
    forward_truth = ref_str[id].strip()[:25]
    if result[i] == -1:
        miss += 1
        continue
    else:
        forward_predice = primers[int(result[i])][0]
        acc.append(forward_predice == forward_truth)
        if forward_predice != forward_truth:
            wrong_distance.append(distance[i])

bins = np.linspace(20, 100, 40)

print('ACC: ', np.mean(acc))
print('DROP: ', miss / len(result))

#plt.hist(wrong_distance, bins=bins, color='red', alpha=0.5, label='wrong', density=True)
#plt.hist(distance, bins=bins, color='blue', alpha=0.5, label='all', density=True)
#plt.show()
