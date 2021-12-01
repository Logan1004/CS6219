import random
import sys
import csv

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

if __name__ == '__main__':
    input_file = sys.argv[1]
    noise, noise_id = read_amplified_strands(input_file)
    output_file = sys.argv[2]
    store_shuffled_strands(output_file, noise, noise_id)

