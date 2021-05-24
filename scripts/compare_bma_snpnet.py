#!/bin/python

import os
import csv
import pandas as pd

complete = pd.read_table('../output/hla_complete.tsv')
with open('../html/HLA_alleles/ukb_to_asterisk_names.csv', mode='r') as infile:
    reader = csv.reader(infile)
    names_dict = {rows[0]:rows[1] for rows in reader}

count = 0
concordance_tot = 0

for gbe_id in os.listdir('../output/snpnet/'):
    if "out" not in gbe_id:
        bma_alleles = list(complete[(complete['GBE_ID'] == gbe_id) & (complete['BMA_posterior_mean'].notna())].sort_values('BMA_posterior_mean', ascending=False)['allelotype'])
        with open('../output/snpnet/' + gbe_id + '/' + gbe_id + '_snpnet.txt') as f:
            snpnet_alleles = f.read().splitlines()
        to_remove = ['PC' + str(i) for i in range(1, 11)] + ['age', 'sex']
        snpnet_alleles = [allele for allele in snpnet_alleles if allele not in to_remove]
        snpnet_alleles = [allele[:-2] for allele in snpnet_alleles]
        snpnet_alleles = [names_dict[allele] for allele in snpnet_alleles]
        if len(bma_alleles) > 0 and len(snpnet_alleles) > 0:
            count += 1
            print(gbe_id)
            num_alleles = min(len(bma_alleles), len(snpnet_alleles), 10)
            concordance = len(set(bma_alleles[:num_alleles]).intersection(set(snpnet_alleles[:num_alleles])))/float(num_alleles)
            concordance_tot += concordance
            print(concordance)
            print(bma_alleles, snpnet_alleles)
print(concordance_tot/count)
