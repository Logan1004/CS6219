import os
import random
import sys
import csv
import numpy as np
import pickle
import shutil
import numba

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


def store_shuffled_strands(output, noise, noise_id, primer_size):
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
    meta['primer size'] = primer_size
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
    forward_vec = np.zeros(grams_size, dtype=bool)
    backward_vec = np.zeros(grams_size, dtype=bool)
    for index, g in enumerate(grams):
        if g in forward:
            forward_vec[index] = 1
        if g in backward:
            backward_vec[index] = 1
    return np.concatenate((forward_vec, backward_vec), axis=0)


def generate_identifier_fast(grams, forward, backward):
    grams_size = len(grams)
    return_vec = np.zeros(2 * grams_size, dtype=bool)
    for index, g in enumerate(grams):
        if g in forward:
            return_vec[index] = 1
        if g in backward:
            return_vec[index + grams_size] = 1
    return return_vec


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

def concatenate_files_v2(input_files, output_file):
    with open(output_file, 'wb') as wfd:
        for tmp_path in input_files:
            with open(tmp_path, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)


def clear_directory(dir_path):
    for f in os.listdir(dir_path):
        file_path = os.path.join(dir_path, f)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


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
        for i in range(len(clusters)):
            f.write("CLUSTER " + str(i) + "\n")
            for s in clusters[i]:
                f.write(s.strip() + "\n")


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


@numba.njit(fastmath=True,cache=True)
def fast_distance_calculate(primer_m, strand_m):
    ret_mat = np.empty((strand_m.shape[0], primer_m.shape[0]))
    for i in range(strand_m.shape[0]):
        for j in range(primer_m.shape[0]):
            acc = 0
            for k in range(strand_m.shape[1]):
                acc += strand_m[i, k] ^ primer_m[j, k]
            ret_mat[i, j] = acc
            #ret_mat[i, j] = np.sum(strand_m[i, :] ^ primer_m[j, :])
    return ret_mat


def test_equal(docu1, docu2, num_of_file):
    for i in range(num_of_file):
        file = 'result-' + str(i) + '.txt'
        path_1 = os.path.join(docu1, file)
        path_2 = os.path.join(docu2, file)
        with open(path_1, 'r') as f: lines_1 = f.readlines()
        with open(path_2, 'r') as f: lines_2 = f.readlines()
        assert lines_1 == lines_2


if __name__ == '__main__':
    # input_file = sys.argv[1]
    # noise, noise_id = read_amplified_strands(input_file)
    # output_file = sys.argv[2]
    # store_shuffled_strands(output_file, noise, noise_id)

    primer_file = sys.argv[1]
    output_file = sys.argv[2]
    primers = read_primers(primer_file)
    for count in [3, 4, 5]:
        mat = []
        grams = generate_n_grams(count)
        for pair in primers:
            identifier = generate_identifier(grams, pair[0], pair[1])
            mat.append(identifier)
        mat = np.array(mat)
        np.save(output_file + '-' + str(count), mat)

    # Command: Python3 Utils.py 'synthesis_data_6219/amplified_strands.txt' 'synthesis_data_6219/UnderlyingClusters.txt'
    # input_file = sys.argv[1]
    # output_file = sys.argv[2]
    # generate_underlying_cluster(input_file, output_file)

    # Command: Python3 Utils.py 'synthesis_data_6219/amplified_strands.txt' 'synthesis_data_6219/raw_data/output.txt'
    # input_file = sys.argv[1]
    # output_file = sys.argv[2]
    # generate_raw_data_from_amplified(input_file, output_file)
