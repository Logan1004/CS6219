import os

import numpy as np

from Utils import *

docu = 'real_data_6219'
n_gram = 5
primer_length = 25

grams = generate_n_grams(n_gram)
primer_id_path = os.path.join(docu, 'Primer', 'primer_gram-' + str(n_gram) + '.npy')
primer_path = os.path.join(docu, 'Primer', 'primer.txt')
primers = read_primers(primer_path)
primer_id_mat = np.load(primer_id_path)
primer_map = {}
for index in range(len(primers)):
    primer_map[tuple(primers[index])] = index

wrong_path = os.path.join(docu, 'wrong_information.txt')
wrong_info = []
with open(wrong_path, 'r') as f:
    wrong_info = f.readlines()

true_id_list = []
wrong_id_list = []
min_edit_id_list = []
for info in wrong_info[1:]:
    info_arr = info.strip().split(',')
    true_id = primer_map[tuple([info_arr[1], info_arr[2]])]
    wrong_id = primer_map[tuple([info_arr[3], info_arr[4]])]
    true_id_list.append(true_id)
    wrong_id_list.append(wrong_id)

    strand = info_arr[-1]
    forward = strand[:primer_length]
    backward = strand[-primer_length:]
    identifier = generate_identifier(grams, forward, backward)

    dist_vec = np.count_nonzero(np.not_equal(primer_id_mat, identifier), axis=1)
    print(dist_vec)
    # primer_distance = np.count_nonzero(np.not_equal(primer_id_mat[true_id], primer_id_mat[wrong_id]))
    # true_distance = np.count_nonzero(np.not_equal(primer_id_mat[true_id], identifier))
    # wrong_distance = np.count_nonzero(np.not_equal(primer_id_mat[wrong_id], identifier))
    #
    # true_forward_edit = edit_distance(info_arr[1], forward, len(info_arr[1]), len(forward))
    # true_backward_edit = edit_distance(info_arr[2], backward, len(info_arr[2]), len(backward))
    # true_edit = [true_forward_edit, true_backward_edit]
    # wrong_forward_edit = edit_distance(info_arr[3], forward, len(info_arr[3]), len(forward))
    # wrong_backward_edit = edit_distance(info_arr[4], backward, len(info_arr[4]), len(backward))
    # wrong_edit = [wrong_forward_edit, wrong_backward_edit]
    forward_edit = []
    backward_edit = []
    potential = np.argsort(dist_vec)[:4]
    for pot_id in potential:
        p = primers[pot_id]
        tmp_forward_edit = edit_distance(p[0], forward, len(p[0]), len(forward))
        forward_edit.append(tmp_forward_edit)
        tmp_backward_edit = edit_distance(p[1], backward, len(p[1]), len(backward))
        backward_edit.append(tmp_backward_edit)

    forward_edit = np.array(forward_edit)
    backward_edit = np.array(backward_edit)
    total_edit = forward_edit + backward_edit
    min_edit_id_list.append(potential[np.argmin(total_edit)])

    print(min_edit_id_list)
    print(total_edit)
    print(true_id_list)
    print(wrong_id_list)
    break

# true_cases = 0
# for i in range(len(min_edit_id_list)):
#     if min_edit_id_list[i] == true_id_list[i]:
#         true_cases += 1
# print(true_cases / len(min_edit_id_list))