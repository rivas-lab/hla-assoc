#!/bin/bash

ml load plink2

plink2 --bfile /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2 --extract <(grep -Fwf <(cut -d':' -f2 sum_table.tsv | grep -v GBE) /oak/stanford/groups/mrivas/ukbb24983/cal/pgen/ukb24983_cal_cALL_v2.bim  | awk -F'\t' '{if ($1 == "6") {print}}' | cut -f2) --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe --remove /oak/stanford/groups/mrivas/ukbb24983/sqc/w24983_20200204.csv --make-bed --out ../output/hla_snp

ml load plink

plink --bfile ../output/hla_snp --recodeAD --out ../output/hla_snp

awk '{for(i=1;i<=NF;i+=2){printf "%s ",$i;} print ""}' ../output/hla_snp.raw | cut -d' ' -f1,4- > ../output/hla_snp_dosage.ped

head -1 ../output/hla_snp_dosage.ped  | tr ' ' '\n' | cut -d'_' -f1 | tr '\n' ' ' > header && echo "" >> header

cat header <(tail -n +2 ../output/hla_snp_dosage.ped) > tmp && mv tmp ../output/hla_snp_dosage.ped

rm ../output/hla_snp.raw
rm header
