# convert_to_ped.py

import numpy as np
import time
import sys

# take in row with one column for each allele, max=2 per entry and convert it
# to a row with two columns for each allele; 'P' means 'present', 'N' means 
# 'not present'
def double_row(row):
    new_row = np.empty((2 * row.size,), dtype=np.unicode)
    for i in range(row.size):
        if row[i] == 0:
            new_row[2*i] = 'N'
            new_row[2*i+1] = 'N'
        elif row[i] == 1:
            new_row[2*i] = 'P'
            new_row[2*i + 1] = 'N'
        elif row[i] == 2:
            new_row[2*i] = 'P' 
            new_row[2*i + 1] = 'P'
        elif row[i] == -1:
            new_row[2*i] = '0'
            new_row[2*i + 1] = '0'
        else:
            raise ValueError('rounded array has value other than 0, 1, 2, -1')
    return new_row

def test_double_row(dosage_values, ped_dosage_values):
    nnz_idx = np.nonzero(dosage_values)
    print('nnz_idx {}'.format(nnz_idx))

    for i in range(nnz_idx[0].size):
#        print('dosage value: {} ped values: {} {}'.format(dosage_values[nnz_idx[0][i], nnz_idx[1][i]], ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2],ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2 + 1]))
        if dosage_values[nnz_idx[0][i], nnz_idx[1][i]] == 1:
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2] == 'P')
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2 + 1] == 'N')
        elif dosage_values[nnz_idx[0][i], nnz_idx[1][i]] == 2:
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2] == 'P')
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2 + 1] == 'P')

# rounds value to nearest integer if within threshold; otherwise returns -1
def round_dosages(val, thresh):
    if val < thresh:
        return 0
    elif 1 - thresh < val and val < 1 + thresh:
        return 1
    elif 2 - thresh < val:
        return 2
    else:
        return -1

def test_round_dosages():
    v_round_dosages = np.vectorize(round_dosages)
    test_a = np.array([[1, 1.4, 0.7],[0.1, 0.05, 2],[0.93, 1.43, 0]])
    print(test_a)
    print(v_round_dosages(test_a, 0.2))

def split_string_underscore(string):
    return string.split('_')[0]

# return length-2 list: first entry is the sum of positive entries of row[idx]
# second is sum of negative entries of row[idx]
def sum_pos_neg(row, idx):
    sub_row = row[idx]
    return [np.sum(sub_row[sub_row > 0]), np.sum(sub_row[sub_row < 0])]

# If a locus has no missing data and its sum != 2, replace all nonzero
# entries for that locus and sample with -1 (missing data)
def check_sum(rounded_dosages, loci_variants, loci_names):
    num_wrong = 0 
    v_split_string_underscore = np.vectorize(split_string_underscore)
    loci_variants_split = v_split_string_underscore(loci_variants)
#    rounded_dosages[0][0] = 1
#    rounded_dosages[5][200] = 1
    for locus in loci_names:
        locus_idx = np.where(loci_variants_split == locus)
        locus_sums = np.apply_along_axis(sum_pos_neg, 1, rounded_dosages, locus_idx)
        rounded_wrong = [i for i in range(rounded_dosages.shape[0]) if locus_sums[i][0] > 2 or (locus_sums[i][1] == 0 and locus_sums[i][0] != 2)]
        num_wrong += len(rounded_wrong)
        for i in rounded_wrong:
             locus_entries = rounded_dosages[i,locus_idx]
             locus_entries[np.nonzero(locus_entries)] = -1
             rounded_dosages[i, locus_idx] = locus_entries
    return rounded_dosages, num_wrong
        
def main():

    t0 = time.time()

    prefix = ''

    if len(sys.argv) == 2:
        prefix = sys.argv[1]

    test = False
    create_rounded = True
    create_ped = False

    num_lines = 1000

    if prefix == '':
    #    dosage_file_name = 'test_' + str(num_lines) + '.txt'
    #    fam_file_name = 'test_' + str(num_lines) + '.fam'
    #    ped_file_name = 'test_' + str(num_lines) + '.ped'
        dosage_file_name = '/scratch/PI/mrivas/ukbb/24983/hla/ukb_hla_v2.txt'
        fam_file_name = '/scratch/PI/mrivas/ukbb/24983/fam/ukb2498_cal_v2_s488374.fam'
        ped_file_name = '/output/convert_to_ped/ukb_hla_v2.ped'       
        prefix = dosage_file_name.split('/')[-1].split('.')[0]

    else:
        dosage_file_name = prefix + '.txt'
        fam_file_name = prefix + '.fam'
        ped_file_name = prefix + '.ped'


    print('begin: {}'.format(time.time() - t0))

    # read in rounded dosage values
    dosage_values = np.loadtxt(dosage_file_name, skiprows = 1)
    print('read dosage: {}'.format(time.time() - t0))
    
    dosage_file = open(dosage_file_name, 'r')
    loci_variants = np.array(dosage_file.readline().split())
    dosage_file.close()

    fam_info = np.loadtxt(fam_file_name, dtype = 'string_', usecols = range(5))

    print('read fam: {}'.format(time.time() - t0))

    loci_names = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5']

    thresh = 0.1
    v_round_dosages = np.vectorize(round_dosages)
    rounded_dosages = v_round_dosages(dosage_values, thresh)
 
    # creates without header
    if create_rounded:
        np.savetxt("output/convert_to_ped/" + prefix + '_rounded.txt', rounded_dosages, delimiter = '\t', fmt = '%s')

    if create_ped:

        print('rounded dosages: {}'.format(time.time() - t0))
    
        rounded_dosages, num_wrong = check_sum(rounded_dosages, loci_variants, loci_names)
    
        print('checked sum: {}'.format(time.time() - t0))
        print('rounded {} entries wrong'.format(num_wrong))
    
        ped_dosage_values = np.apply_along_axis(double_row, 1, rounded_dosages)
    
        print('converted to ped values: {}'.format(time.time() - t0))
    
        zero_col = np.zeros((fam_info.shape[0],), dtype = int)
        all_ped_data = np.transpose(np.vstack((np.transpose(fam_info), zero_col, np.transpose(ped_dosage_values))))
        np.savetxt(ped_file_name, all_ped_data, delimiter = ' ', fmt = '%s')
        print('done: {}'.format(time.time() - t0))
        if test:
            test_double_row(rounded_dosages, ped_dosage_values)
            test_round_dosages()
    
#        check_sum(rounded_dosages, loci_variants, loci_names)
main()
