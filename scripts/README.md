# Script file descriptions

* `reg_test.R`: performs all R regressions as well as BMA
* `write_jobs.py`: submits multiple jobs (for regressions, BMA, etc)
* `write_jobs2.py`: submits multiple jobs (for regressions, BMA, etc). This
  version used to be in the top of the repo (where it was `write_jobs.py`), but
Chris moved it here on 2019/2/2. I'm not sure which of this file or the
`write_jobs.py` in this directory is the most recent version.
* `additivityTest.R`: Chris's example R script
* `convert_to_ped.py`: creates `.ped` file from dosage file for running PLINK
* `count_missing.R`: check which values are missing after rounding dosage file
* `create_fam.py`: creates `.fam` file for running PLINK
* `create_map.py`: creates `.map` file for running PLINK
* `examine_output.R`: puts data from `.rds` object into a text file
* `hap_post_plot.py`: Creates tables of posterior probabilities from BMA output
* `haplotype_names.py`: creates file to translate from literature allelotype
* names to those in our data files
* `make_dosage_rds.R`: Removes unnecessary rows of the dataframe, saves as rds
* (also saves covar version)
* `plot_log.py`: does nothing
* `print_rds.R`: Creates text file from `.rds`
* `process_output.py`: parses R and PLINK regression output to create file with p
values and odds ratios
