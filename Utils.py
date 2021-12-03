import os
import random
import sys
import csv
import numpy as np
import pickle
import shutil

def read_amplified_strands(input_file):
    counter = 0
    noise = []
    noise_id = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            arr = lines[i].strip()
            if arr.isnumeric():
                for j in range(int(arr)):
                    noise.append(lines[i+2+j].strip())
                    noise_id.append(counter)
                counter += 1
    return noise, noise_id


def store_shuffled_strands(output, noise, noise_id):
    index = list(range(len(noise)))
    random.shuffle(index)
    with open(output, 'w') as f:
        for id in index:
            f.write(noise[id] + '\n')

    output2 = output[:-4] + '-index.csv'
    with open(output2, 'w') as f:
        csvwriter = csv.writer(f)
        for id in index:
            csvwriter.writerow([noise_id[id]])

    output3 = output[:-4] + '-meta.pkl'
    meta = {}
    meta['strand size'] = len(noise)
    with open(output3, 'wb') as f:
        pickle.dump(meta, f)


def generate_n_grams(n, alpha=['A', 'T', 'C', 'G']):
    ret_arr = []
    if n == 1:
        ret_arr = alpha
    else:
        tmp_arr = generate_n_grams(n-1)
        for a in alpha:
            ret_arr.extend([a + item for item in tmp_arr])
    return ret_arr


def generate_identifier(grams, forward, backward):
    grams_size = len(grams)
    forward_vec = np.zeros(grams_size)
    backward_vec = np.zeros(grams_size)
    for index, g in enumerate(grams):
        if g in forward:
            forward_vec[index] = 1
        if g in backward:
            backward_vec[index] = 1
    return np.concatenate((forward_vec, backward_vec), axis=0)


def read_primers(primer_file):
    primers = []
    with open(primer_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            arr = line.strip().split(' ')
            primers.append([arr[0], arr[1]])
    return primers


def get_payload_with_primer(strand, primer_pair, primer_length):
    forward = primer_pair[0]
    backward = primer_pair[1]
    forward_pos = [primer_length-1, primer_length, primer_length+1]
    forward_pos_editD = []
    backward_pos = [-primer_length+1, -primer_length, -primer_length-1]
    backward_pos_editD = []

    for i in range(len(forward_pos)):
        forward_substr = strand[: forward_pos[i]]
        forward_pos_editD.append(edit_distance(forward, forward_substr, len(forward), len(forward_substr)))
        backward_substr = strand[backward_pos[i]:]
        backward_pos_editD.append(edit_distance(backward, backward_substr, len(backward), len(backward_substr)))

    f_pos = forward_pos[np.argmin(forward_pos_editD)]
    b_pos = backward_pos[np.argmin(backward_pos_editD)]
    return strand[f_pos: b_pos]


def edit_distance(str1, str2, m, n):
    dp = [[0 for x in range(n + 1)] for x in range(m + 1)]

    # Fill d[][] in bottom up manner
    for i in range(m + 1):
        for j in range(n + 1):

            if i == 0:
                dp[i][j] = j  # Min. operations = j

            elif j == 0:
                dp[i][j] = i  # Min. operations = i

            elif str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]

            else:
                dp[i][j] = 1 + min(dp[i][j - 1],  # Insert
                                   dp[i - 1][j],  # Remove
                                   dp[i - 1][j - 1])  # Replace

    return dp[m][n]

def concatenate_files(input_files, output_file, docu):
    output_path = os.path.join(docu, output_file)
    with open(output_path, 'wb') as wfd:
        for file in input_files:
            tmp_path = os.path.join(docu, file)
            with open(tmp_path, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)


def clear_directory(dir_path):
    for f in os.listdir(dir_path):
        path = os.path.join(dir_path, f)
        os.remove(path)


def generate_underlying_cluster(input_file, output_file):
    clusters = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            arr = lines[i].strip()
            if arr.isnumeric():
                tmp_clu = []
                for j in range(int(arr)):
                    tmp_clu.append(lines[i + 2 + j].strip())
                clusters.append(tmp_clu)

    with open(output_file, 'w') as f:
        csvwriter = csv.writer(f)
        for i in range(len(clusters)):
            csvwriter.writerow(['CLUSTER ' + str(i)])
            for s in clusters[i]:
                csvwriter.writerow([s])


def generate_raw_data_from_amplified(input_file, output_file):
    raw_data = []
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for i in range(len(lines)):
            arr = lines[i].strip()
            if arr.isnumeric():
                raw_data.append([lines[i+1].strip(),int(arr)])

    with open(output_file, 'w') as f:
        csvwriter = csv.writer(f)
        for r in raw_data:
            tmp = r[0] + '_' + str(r[1])
            csvwriter.writerow([tmp])


if __name__ == '__main__':
    # input_file = sys.argv[1]
    # noise, noise_id = read_amplified_strands(input_file)
    # output_file = sys.argv[2]
    # store_shuffled_strands(output_file, noise, noise_id)

    # primer_file = sys.argv[1]
    # output_file = sys.argv[2]
    # primers = read_primers(primer_file)
    # for count in [6, 7]:
    #     mat = []
    #     grams = generate_n_grams(count)
    #     for pair in primers:
    #         identifier = generate_identifier(grams, pair[0], pair[1])
    #         mat.append(identifier)
    #     mat = np.array(mat)
    #     np.save(output_file + '-' + str(count), mat)

    # Command: Python3 Utils.py 'real_data_6219/amplified_strands.txt' 'real_data_6219/UnderlyingClusters.txt'
    # input_file = sys.argv[1]
    # output_file = sys.argv[2]
    # generate_underlying_cluster(input_file, output_file)

    # Command: Python3 Utils.py 'synthesis_data_6219/amplified_strands.txt' 'synthesis_data_6219/raw_data/output.txt'
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    generate_raw_data_from_amplified(input_file, output_file)
