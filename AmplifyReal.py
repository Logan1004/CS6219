import os
import numpy as np
from Utils import *
from Generator import Generator
import matplotlib.pyplot as plt
np.random.seed(6219)
random.seed(6219)


file_name = 'outpout.txt'

data_docu = '0.03-var'
error_rate = 0.03
default_time = 'gamma'
append_primer = False

primer_size = 30
data_path = os.path.join(data_docu, 'raw_data', file_name)

reference = []
times = []
with open(data_path, 'r') as f:
    lines = f.readlines()
    for line in lines:
        arr = line.strip().split('_')
        reference.append(arr[0])
        if len(arr) > 1:
            times.append(int(arr[1]))
        else:
            if default_time == 'gamma':
                times.append(np.maximum(int(np.random.gamma(1.975, 1, 1) * 83.377 - 1.138), 0))
            else:
                times.append(int(default_time))


noise = []
if append_primer:
    forward = 'AAGGCAAGTTGTTACCAGCA'
    backward = 'TGCGACCGTAATCAAACCAA'
else:
    forward = ''
    backward = ''

with open(data_path, 'w') as f:
    for i in range(len(reference)):
        ref_strand = forward + reference[i] + backward
        tmp = []
        for _ in range(times[i]):
            tmp.append(Generator.amplify_strand(ref_strand, error_rate, ['A', 'T', 'C', 'G']))
        noise.append(tmp)
        f.write(ref_strand + '\n')

output_file = 'amplified_strands.txt'
output_path_1 = os.path.join(data_docu, output_file)
counter = []
with open(output_path_1, 'w') as f:
    for i in range(len(reference)):
        num = times[i]
        ref_strand = forward + reference[i] + backward
        f.write(str(num) + '\n')
        f.write(ref_strand + '\n')
        for n in noise[i]:
            f.write(n + '\n')
            counter.append(len(n))
print(np.mean(counter))

# output_path_3 = os.path.join(data_docu, 'UnderlyingClusters.txt')
# with open(output_path_3, 'w') as f:
#     i = 0
#     for cluster in noise:
#         f.write("CLUSTER " + str(i) + "\n")
#         i = i + 1
#         for s in cluster:
#             f.write(s.strip() + "\n")

noise, noise_id = read_amplified_strands(output_path_1)
output_path_2 = os.path.join(data_docu, 'dna_pool.txt')
store_shuffled_strands(output_path_2, noise, noise_id, primer_size)

output_path_3 = os.path.join(data_docu, 'UnderlyingClusters.txt')
generate_underlying_cluster(output_path_1, output_path_3)

output_path_4 = os.path.join(data_docu, 'NoisyStrands.txt')
with open(output_path_4, 'w') as out_file:
    in_file = open(output_path_3, "r")
    l1 = in_file. readlines()
    in_file.close()
    random.shuffle(l1)
    for l in l1:
      if l[0:7] != "CLUSTER":
         out_file.write(l.strip()+"\n")
    out_file.close()

if forward != '':
    for primer_size in [1, 10, 20, 30]:
        subfolder = 'Primer_' + str(primer_size)
        print(subfolder)
        os.makedirs(os.path.join(data_docu, subfolder), exist_ok=True)
        primer_file = os.path.join(data_docu, subfolder, 'primer.txt')
        output_file = os.path.join(data_docu, subfolder, 'primer_gram')
        primers = read_primers(primer_file)
        for count in [3, 4, 5]:
            mat = []
            grams = generate_n_grams(count)
            for pair in primers:
                identifier = generate_identifier(grams, pair[0], pair[1])
                mat.append(identifier)
            mat = np.array(mat)
            np.save(output_file + '-' + str(count), mat)

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
#
# primer_list = set()
# with open(data_path, 'r') as f:
#     lines = f.readlines()
#     for i in range(len(lines)):
#         arr = lines[i].strip().split('_')
#         if arr[0] == '>oligos':
#             primer_pair = tuple([arr[2], arr[3]])
#             primer_list.add(primer_pair)
#
# output_docu = 'real_data_6219/Primer'
# output_file = 'primer.txt'
# output_path = os.path.join(output_docu, output_file)
# with open(output_path, 'w') as f:
#     for pair in primer_list:
#         f.write(pair[0] + ' ' + pair[1] + '\n')





