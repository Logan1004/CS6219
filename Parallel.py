import csv
import time
import pickle
import numpy as np
import itertools
import argparse
import concurrent.futures
from Utils import *
import os


def process_func_fast(docu, start, end, n_gram, primer_length, low_threshold, high_threshold, pid, ext_version, primer_size):
    stat_record = dict()
    stat_record['pid'] = pid
    start_time = time.time()

    input_path = os.path.join(docu, 'dna_pool.txt')
    input_strands = []
    with open(input_path, "r") as text_file:
        for line in itertools.islice(text_file, start, end):
            input_strands.append(line.strip())

    grams = generate_n_grams(n_gram)
    if primer_size == -1:
        primer_folder = 'Primer'
    else:
        primer_folder = 'Primer_' + str(primer_size)
    primer_id_path = os.path.join(docu, primer_folder, 'primer_gram-' + str(n_gram) +'.npy')
    primer_path = os.path.join(docu, primer_folder, 'primer.txt')
    primers = read_primers(primer_path)
    primer_id_mat = np.load(primer_id_path)

    tmp_time = time.time()
    stat_record['initialization'] = tmp_time - start_time
    start_time = tmp_time

    payloads = []
    identifier_mat = np.empty((len(input_strands), primer_id_mat.shape[1]), dtype=bool)
    for strand_id, strand in enumerate(input_strands):
        forward = strand[:primer_length]
        backward = strand[-primer_length:]
        identifier = generate_identifier_fast(grams, forward, backward)
        identifier_mat[strand_id, :] = identifier

    dist_mat = fast_distance_calculate(primer_id_mat, identifier_mat)
    ret_arr = np.argmin(dist_mat, axis=1)

    tmp_time = time.time()
    stat_record['preprocess'] = tmp_time - start_time
    start_time = tmp_time

    for index in range(len(ret_arr)):
        ret_id = ret_arr[index]
        strand = input_strands[index]
        if ret_id == -1 or ext_version == 0:
            payloads.append(strand[primer_length: -primer_length])
        else:
            payloads.append(get_payload_with_primer(strand, primers[ret_id], primer_length))

    tmp_time = time.time()
    stat_record['extraction'] = tmp_time - start_time
    start_time = tmp_time

    output_file_id = os.path.join(docu, 'Output', 'output-' + str(pid) +'-id.csv')
    ret_arr = np.array(ret_arr)
    sorted_idx = np.argsort(ret_arr)
    ret_arr = ret_arr[sorted_idx]
    payloads = np.array(payloads)[sorted_idx]

    with open(output_file_id, 'w') as f:
        for a in ret_arr:
            f.write(str(a) + '\n')

    output_file_id = os.path.join(docu, 'Output', 'output-' + str(pid) +'-payload.txt')
    with open(output_file_id, 'w') as f:
        for s in payloads:
            f.write(s + '\n')

    tmp_time = time.time()
    stat_record['output'] = tmp_time - start_time

    return stat_record

def process_func(docu, start, end, n_gram, primer_length, low_threshold, high_threshold, pid, ext_version, primer_size):
    stat_record = dict()
    stat_record['pid'] = pid
    start_time = time.time()

    input_path = os.path.join(docu, 'dna_pool.txt')
    input_strands = []
    with open(input_path, "r") as text_file:
        for line in itertools.islice(text_file, start, end):
            input_strands.append(line.strip())

    grams = generate_n_grams(n_gram)
    if primer_size == -1:
        primer_folder = 'Primer'
    else:
        primer_folder = 'Primer_' + str(primer_size)
    primer_id_path = os.path.join(docu, primer_folder, 'primer_gram-' + str(n_gram) +'.npy')
    primer_path = os.path.join(docu, primer_folder, 'primer.txt')
    primers = read_primers(primer_path)
    primer_id_mat = np.load(primer_id_path)

    tmp_time =  time.time()
    stat_record['initialization'] = tmp_time - start_time
    start_time = tmp_time

    ret_arr = []
    payloads = []
    for strand_id, strand in enumerate(input_strands):
        forward = strand[:primer_length]
        backward = strand[-primer_length:]
        identifier = generate_identifier(grams, forward, backward)
        dist_vec = np.count_nonzero(np.not_equal(primer_id_mat, identifier), axis=1)
        min_id = np.argmin(dist_vec)
        min_dist = dist_vec[min_id]
        if min_dist < low_threshold:
            ret_arr.append([min_id]) #, min_dist, strand_id+start])
            # payloads.append(get_payload_with_primer(strand, primers[min_id], primer_length))
        elif min_dist > high_threshold:
            ret_arr.append([-1]) # , min_dist, strand_id+start])
            # payloads.append(strand[primer_length, -primer_length])
        else:
            forward_edit = []
            backward_edit = []
            potential = np.argsort(dist_vec)[:3]
            for pot_id in potential:
                p = primers[pot_id]
                tmp_forward_edit = edit_distance(p[0], forward, len(p[0]), len(forward))
                forward_edit.append(tmp_forward_edit)
                tmp_backward_edit = edit_distance(p[1], backward, len(p[1]), len(backward))
                backward_edit.append(tmp_backward_edit)
            total_edit = np.array(forward_edit) + np.array(backward_edit)
            min_id = np.argmin(total_edit)
            ret_arr.append([potential[min_id]]) #, min_dist, strand_id+start])

    tmp_time = time.time()
    stat_record['preprocess'] = tmp_time - start_time
    start_time = tmp_time

    for index in range(len(ret_arr)):
        ret_id = ret_arr[index][0]
        strand = input_strands[index]
        if ret_id == -1 or ext_version == 0:
            payloads.append(strand[primer_length: -primer_length])
        else:
            payloads.append(get_payload_with_primer(strand, primers[ret_id], primer_length))

    tmp_time = time.time()
    stat_record['extraction'] = tmp_time - start_time
    start_time = tmp_time

    output_file_id = os.path.join(docu, 'Output', 'output-' + str(pid) +'-id.csv')
    ret_arr = np.array(ret_arr)
    sorted_idx = np.argsort(ret_arr[:, 0], axis=0)
    ret_arr = ret_arr[sorted_idx]
    payloads = np.array(payloads)[sorted_idx]
    np.savetxt(output_file_id, ret_arr, delimiter=",", fmt='%d')
    # with open(output_file_id, 'w') as f:
    #     writer = csv.writer(f)
    #     writer.writerows(ret_arr)

    output_file_id = os.path.join(docu, 'Output', 'output-' + str(pid) +'-payload.txt')
    np.savetxt(output_file_id, payloads, delimiter=",", fmt='%s')
    # with open(output_file_id, 'w') as f:
    #     for s in payloads:
    #         f.write(s + '\n')

    tmp_time = time.time()
    stat_record['output'] = tmp_time - start_time

    return stat_record

def extract_func(docu, primer_id, process_id, index):
    start_time = time.time()
    input_docu = os.path.join(docu, 'Output')
    output_docu = os.path.join(docu, 'Output_files', str(primer_id))
    os.makedirs(output_docu, exist_ok=True)

    # index = np.genfromtxt(intput_id, delimiter=',', dtype=np.int32)[:, 0]
    # intput_ps_path = os.path.join(docu, 'Output', 'output-' + str(process_id) + '-id.csv')
    # with open(intput_ps_path, 'r') as f:
    #     lines = f.readlines()
    #     index = [int(line.strip()) for line in lines]
    # index = np.array(index, dtype=np.int32)
    position = np.where(index == primer_id)

    if len(position[0]) == 0:
        return time.time() - start_time
    else:
        start_pos = int(np.amin(position))
        end_pos = int(np.amax(position)) + 1

        input_payload = os.path.join(input_docu, 'output-' + str(process_id) + '-payload.txt')
        output_payload = os.path.join(output_docu, 'o-' + str(process_id) + '.txt')
        with open(input_payload, "r") as ext_input_file:
            with open(output_payload, 'w') as ext_output_file:
                ext_output_file.writelines(list(itertools.islice(ext_input_file, start_pos, end_pos)))

        return time.time() - start_time


def reduce_func(docu, primer_id, remove_dir=False):
    start_time = time.time()
    input_docu = os.path.join(docu, 'Output_files', str(primer_id))
    reduce_file_path = os.path.join(docu, 'Output_files', 'result-' + str(primer_id) + '.txt')

    f_list = []
    for f in os.listdir(input_docu):
        file_path = os.path.join(input_docu, f)
        if os.path.isfile(file_path):
            f_list.append(file_path)

    concatenate_files_v2(f_list, reduce_file_path)
    if remove_dir:
        shutil.rmtree(input_docu)

    return time.time() - start_time


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--process', type=int, default=4, help='number of process (default: 4)')
    parser.add_argument('-l', '--primerLen', type=int, default=25, help='primer length (default: 25)')
    parser.add_argument('-w', '--low_threshold', type=int, default=99999, help='low edit distance threshold (default: 99999)')
    parser.add_argument('-v', '--high_threshold', type=int, default=99999, help='high edit distance threshold (default: 99999)')
    parser.add_argument('-g', '--q_grams', type=int, default=5, help='q gram size (default: 5)')
    parser.add_argument('-d', '--docu', type=str, default='real_data_6219', help='input document (default: real_data_6219)')
    parser.add_argument('-e', '--extraction', type=int, default=1, help='0 - naive extraction, 1 - precise extraction (default 1)')
    parser.add_argument('-s', '--primer_size', type=int, default=-1, help='-1 - meta, >1 - customer (default -1)')
    parser.add_argument('-f', '--fast_version', type=int, default=0, help='-1 - fast, 0 - trivial (default 0)')
    args = parser.parse_args()

    data_docu = args.docu
    process_num = int(args.process)
    q_grams = int(args.q_grams)
    primer_length = int(args.primerLen)
    low_threshold = int(args.low_threshold)
    high_threshold = int(args.high_threshold)
    extraction_version = int(args.extraction)
    if args.fast_version == 1:
        use_func = process_func_fast
    else:
        use_func = process_func

    os.makedirs(os.path.join(data_docu, 'Output'), exist_ok=True)
    os.makedirs(os.path.join(data_docu, 'Output_files'), exist_ok=True)
    clear_directory(os.path.join(data_docu, 'Output'))
    clear_directory(os.path.join(data_docu, 'Output_files'))

    meta_data_file = 'dna_pool-meta.pkl'
    meta_data_path = os.path.join(data_docu, meta_data_file)
    with open(meta_data_path, "rb") as f: meta = pickle.load(f)
    if int(args.primer_size) != -1:
        meta['primer size'] = int(args.primer_size)

    cpu_count = min(os.cpu_count(), process_num)
    strands_per_cpu = int(np.ceil( meta['strand size'] / cpu_count))

    start = time.time()
    future_list = []
    with concurrent.futures.ProcessPoolExecutor(cpu_count) as executor:
        for pid in range(cpu_count):
            begin = pid * strands_per_cpu
            end = min(begin + strands_per_cpu, meta['strand size'])
            tmp_fut = executor.submit(use_func, data_docu, begin, end, q_grams, primer_length,
                                      low_threshold, high_threshold, pid, extraction_version, int(args.primer_size))
            future_list.append(tmp_fut)
    # for f in future_list:
    #     print(f.result())
    parallel_time = time.time() - start
    # print('Parallel Processing Time:', parallel_time)

    t_future_list_ext = []
    t_future_list_red = []
    start = time.time()
    with concurrent.futures.ProcessPoolExecutor(cpu_count) as ext_executor:
        for ps in range(cpu_count):
            intput_ps_path = os.path.join(data_docu, 'Output', 'output-' + str(ps) + '-id.csv')
            with open(intput_ps_path, 'r') as f:
                lines = f.readlines()
                index = [int(line.strip()) for line in lines]
            index = np.array(index, dtype=np.int32)
            for fid in range(meta['primer size']):
                tmp_fut = ext_executor.submit(extract_func, data_docu, fid, ps, np.copy(index))
                t_future_list_ext.append(tmp_fut)
    check_point1 = time.time() - start

    with concurrent.futures.ProcessPoolExecutor(cpu_count) as red_executor:
        for fid in range(meta['primer size']):
            tmp_fut = red_executor.submit(reduce_func, data_docu, fid, True)
            t_future_list_red.append(tmp_fut)
    # for fid in range(meta['primer size']):
    #     reduce_func(data_docu, fid, True)
    reduce_time = time.time() - start
    # print(reduce_time)

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
        writer.writerow(['Total', -1, parallel_time, -1, reduce_time])
        for fut in future_list:
            tmp_stat = fut.result()
            writer.writerow([tmp_stat['pid'], tmp_stat['initialization'],
                             tmp_stat['preprocess'], tmp_stat['extraction'],
                             tmp_stat['output'] ])

    print("----- Total time: ", parallel_time + reduce_time, " -----")
    #
    ext_time_per = []
    for f in t_future_list_ext:
        ext_time_per.append(f.result())
    ext_time_path = os.path.join(data_docu, 'extract.pickle')
    with open(ext_time_path, 'wb') as f:
       pickle.dump(ext_time_per, f)
    red_time_per = []
    for f in t_future_list_red:
        red_time_per.append(f.result())
    # print('EXT TIME: ', ext_time_per)
    # print('RED TIME: ', red_time_per)
    # print('MAX EXT: ', np.amax(ext_time_per))
    # print('EXTRACTION TIME: ', check_point1)