import os
import random
import csv
import argparse
from Utils import *
import numpy as np
np.random.seed(6219)
file_size_range = [700, 1400]

class Generator:
    def __init__(self, num_of_file, error_rate, data_docu, duplicate):
        self.data_document = os.path.join(data_docu, 'Primer')
        self.primer_file = 'primer.txt'
        self.alpha = ['A', 'T', 'C', 'G']
        self.num_of_file = num_of_file
        self.err_rate = error_rate
        if duplicate == 'gamma':
            self.gamma = True
        else:
            self.gamma = False
            self.duplicate = int(duplicate)

        self.primers_arr = []
        self.read_primer() # read primers
        self.original_files = []
        self.generate_files() # generate files
        self.amplified_pool = []
        self.amplified_strand_info = []
        self.info_map_strand = {}
        self.amplify_dna_pool() # amplify the dna pool
        # self.shuffled_index = []
        # self.shuffled_pool = []
        # self.shuffle_dna_pool()  # amplify the dna pool

    def read_primer(self):
        file = os.path.join(self.data_document, self.primer_file)
        with open(file) as f: lines = f.readlines()
        for line in lines:
            line = line.strip().split(' ')
            self.primers_arr.append(tuple(line))

    def generate_files(self):
        for i in range(self.num_of_file):
            self.original_files.append(self.generate_single_file(i))

    def get_index(self, index, base=6):
        # 6 bases length index
        index_base = ['A' for _ in range(base)]
        for i in range(base):
            index_base[base-i-1] = self.alpha[index % 4]
            index = index // 4
        return ''.join(index_base)

    def generate_single_file(self, primer_id, payload_size=100):
        #file_size = random.randint(file_size_range[0], file_size_range[1])
        file_size = 2000 #65534
        ret_arr = []
        forward_primer = self.primers_arr[primer_id][0]
        backward_primer = self.primers_arr[primer_id][1]
        for index in range(file_size):
            tmp = forward_primer + self.get_index(index)
            for _ in range(payload_size):
                tmp += random.choice(self.alpha)
            tmp += backward_primer
            ret_arr.append(tmp)
        return ret_arr

    def amplify_dna_pool(self):
        for f_id, file in enumerate(self.original_files):
            for s_id, strand in enumerate(file):
                info = tuple([f_id, s_id])
                self.info_map_strand[info] = []

                if self.gamma:
                    duplicate = np.maximum(int(np.random.gamma(1.975, 1, 1) * 83.377 - 1.138), 0)
                else:
                    duplicate = self.duplicate

                for _ in range(duplicate):
                    tmp = Generator.amplify_strand(strand, self.err_rate, self.alpha)
                    self.amplified_pool.append(tmp)
                    self.amplified_strand_info.append(info)
                    self.info_map_strand[info].append(len(self.amplified_pool)-1)

    @staticmethod
    def amplify_strand(strand, err, alpha):
        tmp = ''
        for b in strand:
            tie = random.random()
            if tie < err / 3:
                # insert before
                tmp += random.choice(alpha) + b
            elif tie < 2 * err / 3:
                # delete
                tmp += ''
            elif tie < err:
                # substitute
                tmp += random.choice(alpha)
            else:
                tmp += b
        return tmp

    def shuffle_dna_pool(self):
        pool_size = len(self.amplified_pool)
        self.shuffled_index = list(range(pool_size))
        self.shuffled_pool = ['' for _ in range(pool_size)]
        random.shuffle(self.shuffled_index)
        for i in range(pool_size):
            s_id = self.shuffled_index[i]
            self.shuffled_pool[i] = self.amplified_pool[s_id]

    def store_shuffled_pool(self, output_file):
        pool_size = len(self.amplified_pool)
        with open(output_file, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['strand', 'file id', 'pos in file', 'ground truth'])
            for i in range(pool_size):
                strand = self.shuffled_pool[i]
                original_id = self.shuffled_index[i]
                info = self.amplified_strand_info[original_id]
                ground_truth = self.original_files[info[0]][info[1]]
                writer.writerow([strand, info[0], info[1], ground_truth])

    def store_clustered_pool(self, output_file):
        with open(output_file, 'w') as f:
            for key, value in self.info_map_strand.items():
                ref_strand = self.original_files[key[0]][key[1]]
                strand_size = len(value)
                f.write(str(strand_size) + '\n')
                f.write(ref_strand + '\n')
                for index in value:
                    f.write(self.amplified_pool[index]+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--num_of_files', type=int, default=1, help='Number of files generated, ranging from 0 to 35 (default: 1)')
    parser.add_argument('-e', '--error_rate', type=float, default=0.0, help='Error rate, ranging from 0 to 1 (default: 0.0)')
    parser.add_argument('-d', '--data_document', type=str, default='synthesis_data_6219', help='input data document (default: synthesis_data_6219)')
    parser.add_argument('-p', '--duplicate', type=str, default='gamma', help='data duplicate (default: gamma)')
    args = parser.parse_args()

    data_docu = args.data_document

    output_path_1 = os.path.join(data_docu, 'amplified_strands.txt')
    gen = Generator(num_of_file=args.num_of_files, error_rate=args.error_rate,
                    data_docu=data_docu, duplicate=args.duplicate)
    gen.store_clustered_pool(output_path_1)

    noise, noise_id = read_amplified_strands(output_path_1)
    output_path_2 = os.path.join(data_docu, 'dna_pool.txt')
    store_shuffled_strands(output_path_2, noise, noise_id, 35)

    output_path_3 = os.path.join(data_docu, 'UnderlyingClusters.txt')
    generate_underlying_cluster(output_path_1, output_path_3)

    os.makedirs(os.path.join(data_docu, 'raw_data'), exist_ok=True)
    output_path_4 = os.path.join(data_docu, 'raw_data', 'output.txt')
    generate_raw_data_from_amplified(output_path_1, output_path_4)


