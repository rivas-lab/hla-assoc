# hla-assoc

## Cloning

Clone using `git lfs clone git@github.com:rivas-lab/hla-assoc.git`.
If you don't use LFS to clone, you may need to to do `git lfs pull` after you
clone to pull in the files tracked by LFS.

## Git LFS

The files in `output` directory are tracked with git LFS.

## Directories

### `output`

This directory is tracked by Git LFS and contains one output directory per
notebook or script from the `notebooks` and `scripts` directories. These
output files are small enough to reasonably track with Git LFS and should not
contain identifiable information or other sensitive information and should be
safe to share publicly, though care should be taken to double-check that before
sharing.

### `private_data`

This directory is **not** tracked by Git and contains one output directory per
notebook or script from the `notebooks` and `scripts` directories. These files
may contain sensitive information, be too large to track with Git LFS, or 
otherwise don't need to go into the output directory at this time.

## Notes

### Rounding for dosage in .ped file

Round every dosage within some threshold to the nearest integer, and replace
dosages outside of the threshold with the "missing data" tag. For example, if
thresh = 0.1, then any value in (0, 0.1) would be rounded to 0, any value in
(0.9, 1.1) would be rounded to 1, and any value in (1.9, 2) would be rounded to 
two. Any values outside these intervals are missing data (-1).

We expect each locus to have entries that sum to 2 for each sample. If the
sum of the entries is >2, or if the sum doesn't equal 2 and yet there is no
missing data, replace all nonzero entries for that sample/locus combo with the
missing data symbol.

I have been using thresh = 0.1 so far.

### PLINK vs R regressions

We ran PLINK and R both with the white British subset of samples, both with 
age, sex, and PC1-PC4 as covariates. For the PLINK run we used a `.bed` file
based on the rounded dosages. (Chris/Guhan: PLINK results on oak at `/oak/stanford/groups/mrivas/users/jolivier/repos/hla-assoc/data/PLINK_results`.)

We ran three regressions in R on each of the haplotypes with >5 nonzero entries
(312 out of 362): an additive logistic regression on the dosages,
additive logistic regression on the integer dosages (rounded as described
above), and a "Dosage" logistic regression where we convert the integer dosage
values to factors. We only perform this regression on haplotypes with >2
distinct dosage values (151 haplotypes).

Here is the comparison of the log odds ratio from the plink run and the
additive regression on integer dosages in R (not including the plink runs that
used firth):

![alt text](https://github.com/rivas-lab/hla-assoc/blob/master/plots/out_dosage_add_plot.png)

We can see that there are a few outliers.

### Testing concordance for Celiac disease

We identified DQA1\_501 and DQB1\_201 as haplotypes related to Celiac disease
from the literature. Here are the plots of the sorted -log<sub>10</sub> p
values and the log odds ratios for each haplotype. The haplotypes of interest
are highlighted in yellow. We can see that their position in these plots makes
sense.

#### Odds ratio

##### PLINK

![alt text](https://github.com/rivas-lab/hla-assoc/blob/master/plots/OR_PLINK_plot.png)

##### R

![alt text](https://github.com/rivas-lab/hla-assoc/blob/master/plots/OR_R_plot.png)


#### P value

The p values from the R regressions are cut off at 2e-16.

##### PLINK

![alt text](https://github.com/rivas-lab/hla-assoc/blob/master/plots/log10_pval_PLINK_plot.png)

##### R
![alt text](https://github.com/rivas-lab/hla-assoc/blob/master/plots/log10_pval_R_plot.png)

