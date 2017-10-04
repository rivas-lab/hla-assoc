# create_map.py

# first column, chromosome, is '27' for all because we don't have chromosome number
# second column is the varient identifier from header of original dosage file
# third column, base-pair coordinate, is arbitrarily labeled with the index

import numpy as np

def main():
    dosage_file_name = '/scratch/PI/mrivas/ukbb/24983/hla/ukb_hla_v2.txt'
    map_file_name = 'ukb_hla_v2.map'

    dosage_file = open(dosage_file_name, 'r')
    loci_variants = np.array(dosage_file.readline().split())
    dosage_file.close()

    col0 = np.full(loci_variants.shape, 27)
    col2 = np.arange(loci_variants.size)

    map_data = np.transpose(np.vstack((col0, loci_variants, col2)))

    np.savetxt(map_file_name, map_data, delimiter = ' ', fmt = '%s')

main()
