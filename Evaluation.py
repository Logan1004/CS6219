import csv

from Utils import *
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--docu', type=str, default='real_data_6219', help='input document (default: real_data_6219)')
args = parser.parse_args()

docu = args.docu
primer_path = os.path.join(docu, 'Primer', 'primer.txt')
primers = read_primers(primer_path)
primer_length = len(primers[0][0])
reference = os.path.join(docu, 'raw_data', 'outpout.txt')
with open(reference, 'r') as f: ref_str = f.readlines()

result_path = os.path.join(docu, 'Output', 'output-id.csv')
result = []
lines = None
distance = []
strand_id = []
wrong_distance = []
wrong_case = []
wrong_id = []
tricky_case = []
true_distance = []
with open(result_path, 'r') as f:
    lines = f.readlines()
    for line in lines:
        arr = line.strip().split(',')
        result.append(int(arr[0]))
        distance.append(int(arr[1]))
        strand_id.append(int(arr[-1]))

index_path = os.path.join(docu, 'dna_pool-index.csv')
with open(index_path, 'r') as f: index = f.readlines()
index = [int(a.strip()) for a in index]

acc = []
miss = 0
for i in range(len(result)):
    s_id = strand_id[i]
    id = index[s_id]
    forward_truth = ref_str[id].strip()[:primer_length]
    if result[i] == -1:
        miss += 1
        continue
    else:
        forward_predice = primers[int(result[i])][0]
        acc.append(forward_predice == forward_truth)
        if forward_predice != forward_truth:
            wrong_distance.append(distance[i])
            wrong_case.append(lines[i])
            wrong_id.append(i)
        else:
            true_distance.append(distance[i])
            if distance[i] >60:
                tricky_case.append(lines[i])

bins = np.linspace(20, 100, 40)

print('ACC: ', np.mean(acc))
print('DROP: ', miss / len(result))

# plt.hist(wrong_distance, bins=bins, color='red', alpha=0.5, label='wrong', density=True)
# plt.hist(distance, bins=bins, color='blue', alpha=0.5, label='true', density=True)
# plt.show()

# for c in wrong_case[:10]:
#     print(c)
# print('################')
# for c in tricky_case[:10]:
#     print(c)

pool_path = os.path.join(docu, 'dna_pool.txt')
with open(pool_path, 'r') as f:
    pool = f.readlines()
primer_path = os.path.join(docu, 'Primer', 'primer.txt')
primers = read_primers(primer_path)
wrong_info_path = os.path.join(docu, 'wrong_information.txt')
with open(wrong_info_path, 'w') as f:
    csvwriter = csv.writer(f)
    csvwriter.writerow(['index', 'truth forward', 'truth backward', 'predict forward', 'predict backward', 'strand'])
    for i in wrong_id:
        label_id = index[i]
        predict_id = result[i]
        strand = pool[i]
        ref_strand = ref_str[label_id].strip().split('_')[0]
        csvwriter.writerow([i, ref_strand[:primer_length], ref_strand[-primer_length:],
                            primers[predict_id][0], primers[predict_id][1], strand.strip()])