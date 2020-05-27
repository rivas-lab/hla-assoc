#!/bin/python

import numpy as np
import pandas as pd
from statsmodels.stats import multitest
import sys

def BY_adjust_p(in_file, col, newcol):
    mask = np.isfinite(in_file[col])
    pval_corrected = np.full(in_file[col].shape, np.nan)
    pval_corrected[mask] = multitest.multipletests(in_file[col][mask], is_sorted=False, method='fdr_by')[1]
    in_file[newcol] = pval_corrected
    return in_file

if __name__ == "__main__":
    in_file = pd.read_table(sys.argv[1])
    print(sys.argv[1])
    col = sys.argv[2]
    print(col)
    newcol = sys.argv[3]
    all_p = BY_adjust_p(in_file, col, newcol)
    all_p.to_csv(sys.argv[1].split('.tsv')[0] + '_adj.tsv', sep='\t', index=False)
