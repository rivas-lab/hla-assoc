# convert_to_ped.py

import numpy as np

# take in row with one column for each allele, max=2 per entry and convert it
# to a row with two columns for each allele, max=1 per entry
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
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2] == b'P')
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2 + 1] == b'N')
        elif dosage_values[nnz_idx[0][i], nnz_idx[1][i]] == 2:
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2] == b'P')
            assert(ped_dosage_values[nnz_idx[0][i], nnz_idx[1][i]*2 + 1] == b'P')

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

def sum_pos_neg(row, idx):
    sub_row = row[idx]
    
    #print(sub_row)
    #print(sub_row > 0)
    return [np.sum(sub_row[sub_row > 0]), np.sum(sub_row[sub_row < 0])]

def check_sum(rounded_dosages, loci_variants, loci_names):
    
    v_split_string_underscore = np.vectorize(split_string_underscore)
    loci_variants_split = v_split_string_underscore(loci_variants)
    for locus in loci_names:
        locus_idx = np.where(loci_variants_split == locus)
        print(locus)
    #    print(rounded_dosages[0,loci_idx])
        locus_sums = np.apply_along_axis(sum_pos_neg, 1, rounded_dosages, locus_idx)
        
def main():

    test = False

    # read in rounded dosage values
    dosage_values = np.loadtxt('test_100.txt', skiprows = 1)
    
    dosage_file = open('test_100.txt', 'r')
    loci_variants = np.array(dosage_file.readline().split())
    dosage_file.close()

    fam_info = np.loadtxt('test_100.fam', dtype = 'string_', usecols = range(5))
    loci_names = ['A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1', 'DRB3', 'DRB4', 'DRB5']

    thresh = 0.1
    v_round_dosages = np.vectorize(round_dosages)
    rounded_dosages = v_round_dosages(dosage_values, thresh)
    print(rounded_dosages)
    # don't all sum to 22
#    for i in range(dosage_values.shape[0]):
#        #print(dosage_values[i, np.nonzero(dosage_values[i,:])])
#        print()
#        print(dosage_values[i,:])
#        print('doubled:')
#        print(double_row(dosage_values[i,:]))
##        print(np.sum(dosage_values[i,:]))
##        print(dosage_values[i,:].shape)
    print(dosage_values)


    ped_dosage_values = np.apply_along_axis(double_row, 1, rounded_dosages)
    print(ped_dosage_values)

    print(dosage_values.shape)
    print(ped_dosage_values.shape)

    print(fam_info.shape)
    print(ped_dosage_values.shape)
    zero_col = np.zeros((fam_info.shape[0],), dtype = int)
    print(zero_col.shape)
    all_ped_data = np.transpose(np.vstack((np.transpose(fam_info), zero_col, np.transpose(ped_dosage_values))))
    np.savetxt('test_100.ped', all_ped_data, delimiter = ' ', fmt = '%s')
    if test:
        test_double_row(rounded_dosages, ped_dosage_values)
        test_round_dosages()
    
        check_sum(rounded_dosages, loci_variants, loci_names)
main()
