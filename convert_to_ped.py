# convert_to_ped.py

import numpy as np

# take in row with one column for each allele, max=2 per entry and convert it
# to a row with two columns for each allele, max=1 per entry
def double_row(row):
    new_row = np.empty((2 * row.size,), dtype='string_')
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


def main():

    # read in rounded dosage values
    dosage_values = np.around(np.loadtxt('test_100.txt', skiprows = 1))

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


    ped_dosage_values = np.apply_along_axis(double_row, 1, dosage_values)
    print(ped_dosage_values)

    print(dosage_values.shape)
    print(ped_dosage_values.shape)
    test_double_row(dosage_values, ped_dosage_values)
main()
