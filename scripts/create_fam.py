import numpy as np
import pandas as pd

def decode_strings(col):
    return np.array([x.decode() for x in col])


def main():
    fam_file = '/scratch/PI/mrivas/ukbb/24983/fam/ukb2498_cal_v2_s488374.fam'
    fam_data = np.loadtxt(fam_file, dtype='string_',usecols=range(5))
    zero_col = np.zeros((fam_data.shape[0],), dtype='int')
    full_data = np.transpose(np.vstack((np.transpose(fam_data), zero_col)))
    print(full_data)
    #full_data.astype('<U100')
#    full_data_decode = np.array([x.decode() for x in full_data])
    full_data = np.apply_along_axis(decode_strings, 0, full_data)
    print(full_data)
    np.savetxt('output/create_fam/ukb_hla_v2.fam', full_data, delimiter = ' ', fmt = '%s')


#    fam_df = pd.read_csv(fam_file, delimiter = ' ', usecols=range(5))
#    print(fam_df)
#    fam_df['last']=0
#    print(fam_df)
#    fam_df.to_csv('ukb_hla_v2.fam', sep = ' ', header=None, index=False)
main()
