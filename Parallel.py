import csv
import time
import pickle
import numpy as np
import itertools
import concurrent.futures
from Utils import *
import os


def process_func(docu, start, end, n_gram, primer_length, threshold, pid):
    stat_record = dict()
    stat_record['pid'] = pid
    start_time = time.time()

    input_path = os.path.join(docu, 'dna_pool.txt')
    input_strands = []
    with open(input_path, "r") as text_file:
        for line in itertools.islice(text_file, start, end):
            input_strands.append(line.strip())

    grams = generate_n_grams(n_gram)
    primer_id_path = os.path.join(docu, 'Primer', 'primer_gram-' + str(n_gram) +'.npy')
    primer_path = os.path.join(docu, 'Primer', 'primer.txt')
    primers = read_primers(primer_path)
    primer_id_mat = np.load(primer_id_path)

    tmp_time =  time.time()
    stat_record['initialization'] = tmp_time - start_time
    start_time = tmp_time

    ret_arr = []
    payloads = []
    for strand in input_strands:
        forward = strand[:primer_length]
        backward = strand[-primer_length:]
        identifier = generate_identifier(grams, forward, backward)
        dist_vec = np.count_nonzero(np.not_equal(primer_id_mat, identifier), axis=1)
        min_id = np.argmin(dist_vec)
        min_dist = dist_vec[min_id]
        if min_dist > threshold:
            ret_arr.append([-1, min_dist])
            # payloads.append(strand[primer_length, -primer_length])
        else:
            ret_arr.append([min_id, min_dist])
            # payloads.append(get_payload_with_primer(strand, primers[min_id], primer_length))

    tmp_time = time.time()
    stat_record['preprocess'] = tmp_time - start_time
    start_time = tmp_time

    for index in range(len(ret_arr)):
        ret_id = ret_arr[index][0]
        strand = input_strands[index]
        if ret_id == -1:
            payloads.append(strand[primer_length: -primer_length])
        else:
            payloads.append(get_payload_with_primer(strand, primers[ret_id], primer_length))

    tmp_time = time.time()
    stat_record['extraction'] = tmp_time - start_time
    start_time = tmp_time

    output_file_id = os.path.join(docu, 'Output', 'output-' + str(pid) +'-id.csv')
    with open(output_file_id, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(ret_arr)

    output_file_id = os.path.join(docu, 'Output', 'output-' + str(pid) +'-payload.txt')
    with open(output_file_id, 'w') as f:
        for s in payloads:
            f.write(s + '\n')

    tmp_time = time.time()
    stat_record['output'] = tmp_time - start_time

    return stat_record


if __name__ == '__main__':
    data_docu = 'real_data_6219'
    process_num = 4
    q_grams = 5
    primer_length = 25
    threshold = 70

    clear_directory(os.path.join(data_docu, 'Output'))

    meta_data_file = 'dna_pool-meta.pkl'
    meta_data_path = os.path.join(data_docu, meta_data_file)
    with open(meta_data_path, "rb") as f: meta = pickle.load(f)

    # meta['strand size'] = 800

    cpu_count = min(os.cpu_count(), process_num)
    strands_per_cpu = int(np.ceil( meta['strand size'] / cpu_count))

    start = time.time()
    future_list = []
    with concurrent.futures.ProcessPoolExecutor(cpu_count) as executor:
        for pid in range(cpu_count):
            begin = pid * strands_per_cpu
            end = min(begin + strands_per_cpu, meta['strand size'])
            tmp_fut = executor.submit(process_func, data_docu, begin, end, q_grams, primer_length, threshold, pid)
            future_list.append(tmp_fut)
    parallel_time = time.time() - start
    # print('Parallel Processing Time:', parallel_time)

    start = time.time()
    id_sub_file = ['output-' + str(index) + '-id.csv' for index in range(cpu_count)]
    id_output_file = 'output-id.csv'
    concatenate_files(id_sub_file, id_output_file, os.path.join(data_docu, 'Output'))
    payload_sub_file = ['output-' + str(index) + '-payload.txt' for index in range(2)]
    payload_output_file = 'output-payload.csv'
    concatenate_files(payload_sub_file, payload_output_file, os.path.join(data_docu, 'Output'))
    merge_time = time.time() - start
    # print('Merging Time: ', merge_time)

    # individual processing time
    output_file_stat = os.path.join(data_docu, 'Output', 'running_statistics.csv')
    with open(output_file_stat, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['PID', 'initialization', 'preprocess', 'extraction', 'output'])
        writer.writerow(['Total', -1, parallel_time, -1, merge_time])
        for fut in future_list:
            tmp_stat = fut.result()
            writer.writerow([tmp_stat['pid'], tmp_stat['initialization'],
                             tmp_stat['preprocess'], tmp_stat['extraction'],
                             tmp_stat['output'] ])