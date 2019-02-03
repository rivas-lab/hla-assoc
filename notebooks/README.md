# Notebook Descriptions

* `Check_bim`: nothing important
* `Full_HLA_frequency`: Find frequency of each allelotype
* `check_firth`: identified which regressions were using firth (what the allele
  variant frequencies were in these cases, what the phenotype frequencies were
in these cases)
* `compare_cutoff_pvals`: found number of associations for various allele
  frequency/phenotype frequency cutoffs. Heatmaps of adjusted p values.
Determines which associations will be used for BMA. Makes frequency donut
charts.
* `count_missing`: counts how many of each value (-1,0,1,2) appear in the
  rounded file
* `dosage_exploration`: Histograms of dosage values and number of entries
  outside the threshold (used to choose the cutoff of 0.1 for rounding)
* `frequency_stratification`: Creates tables of posterior probabilities broken
  down by the frequency of the allelotype in the population
* `phe_hap_table`: initial plots of posterior probabilites; switched to plotly
  later
* `plink_vs_R_analysis`: compares the p values and odds ratios from the plink
  and R analyses (with plots)

* `Process Plink Results`: Read raw plink output files and extract results for
  HLA alleles. 
