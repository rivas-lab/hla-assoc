# Files

* `AIC_add.csv`: AIC calculations for the additive allelotype-phenotype models.
Row is phenotype, column is allelotype.

* `AIC_diff.csv`: Difference in AIC calculations between the additive and genotype
allelotype-phenotype models. Row is phenotype, column is allelotype.

* `AIC_gen.csv`: AIC calculations for the genotype allelotype-phenotype models.
Row is phenotype, column is allelotype.

* `additive_assoc_adj_p_all.csv`: All BY-adjusted p-values from additive model single 
variant association testing. Row is phenotype, column is allelotype.

* `additive_assoc_adj_p_sig.csv`: BY-adjusted p-values from additive model single 
variant association testing. Row is phenotype, column is allelotype. Only selected 
rows and columns that had significant associations (p<=0.05).

* `additive_logOR.csv`: All LogORs from additive model single 
variant association testing. Row is phenotype, column is allelotype.

* `additive_se.csv`: All SEs from additive model single 
variant association testing. Row is phenotype, column is allelotype.

* `bma_post_prob_all.csv`: All results from the BMA analysis, which looked at
a subset of 48 phenotypes and 62 allelotypes. Criteria for selection: at least
2 allelotypes significantly associated with phenotype; took only the top 10
most-significant allelotypes per phenotype if phenotype had more than 10 sig.
associations.

* `bma_post_prob_sig_literature.csv`: Information on whether the 60 putative
  causal alleles identified by BMA have been reported in the literature. Guhan
made this table recently and checked all 60.

* `complete_p_logor_se_aic.tsv`: Contains allelotype, phenotype, and [logOR, p, se]
for [additive, genotype_1, genotype_2] models. Also contains AIC for additive vs
genotype models and the difference in AIC.

* `complete_p_logor_se_aic_desc.tsv`: Everything contained in the above with the
addition of annotated human-readable phenotypes.

* `genotype_assoc_adj_p_all_1.csv`: All BY-adjusted p-values from genotype model 
single variant association testing if the person has 1 of the allelotypes. Having either
one or two of the given allelotype are being treated as separate events that aren't 
assumed to be linearly related to each other in this model. Each have their own odds 
ratios, p-values, and SEs. Row is phenotype, column is allelotype.

* `genotype_assoc_adj_p_all_2.csv`: All BY-adjusted p-values from genotype model 
single variant association testing if the person has 2 of the allelotypes. Having either
one or two of the given allelotype are being treated as separate events that aren't 
assumed to be linearly related to each other in this model. Each have their own odds 
ratios, p-values, and SEs. Row is phenotype, column is allelotype.

* `genotype_assoc_adj_p_sig_1.csv`: BY-adjusted p-values from genotype model 
single variant association testing if the person has 1 of the allelotypes. Having either
one or two of the given allelotype are being treated as separate events that aren't 
assumed to be linearly related to each other in this model. Each have their own odds 
ratios, p-values, and SEs. Row is phenotype, column is allelotype. Only selected rows 
and columns that had significant associations (p<=0.05).

* `genotype_assoc_adj_p_sig_2.csv`: BY-adjusted p-values from genotype model 
single variant association testing if the person has 2 of the allelotypes. Having either
one or two of the given allelotype are being treated as separate events that aren't 
assumed to be linearly related to each other in this model. Each have their own odds 
ratios, p-values, and SEs. Row is phenotype, column is allelotype. Only selected rows 
and columns that had significant associations (p<=0.05).

* `genotype_assoc_unadj_p_all_1.csv`: All unadjusted p-values from genotype model 
single variant association testing if the person has 1 of the allelotypes. Having either
one or two of the given allelotype are being treated as separate events that aren't 
assumed to be linearly related to each other in this model. Each have their own odds 
ratios, p-values, and SEs. Row is phenotype, column is allelotype.

* `genotype_assoc_unadj_p_all_2.csv`: All unadjusted p-values from genotype model 
single variant association testing if the person has 2 of the allelotypes. Having either
one or two of the given allelotype are being treated as separate events that aren't 
assumed to be linearly related to each other in this model. Each have their own odds 
ratios, p-values, and SEs. Row is phenotype, column is allelotype.

* `genotype_logOR_1.csv`: All LogORs from genotype model single variant association
testing if the person has 1 of the allelotypes. Having either one or two of the given
allelotype are being treated as separate events that aren't assumed to be linearly 
related to each other in this model. Each have their own odds ratios, p-values, and SEs.
Row is phenotype, column is allelotype.

* `genotype_logOR_2.csv`: All LogORs from genotype model single variant association
testing if the person has 2 of the allelotypes. Having either one or two of the given
allelotype are being treated as separate events that aren't assumed to be linearly 
related to each other in this model. Each have their own odds ratios, p-values, and SEs.
Row is phenotype, column is allelotype.

* `genotype_se_1.csv`: All SEs from genotype model single variant association
testing if the person has 1 of the allelotypes. Having either one or two of the given
allelotype are being treated as separate events that aren't assumed to be linearly 
related to each other in this model. Each have their own odds ratios, p-values, and SEs.
Row is phenotype, column is allelotype.

* `genotype_se_2.csv`: All SEs from genotype model single variant association
testing if the person has 2 of the allelotypes. Having either one or two of the given
allelotype are being treated as separate events that aren't assumed to be linearly 
related to each other in this model. Each have their own odds ratios, p-values, and SEs.
Row is phenotype, column is allelotype.

* `interaction_adj_p_all.csv`: All BY-adjusted p-values from interaction model 
association testing. Row is phenotype, column is pairs of allelotypes shown to be
significantly associated with phenotype (if the cell is full).

* `interaction_adj_p_sig.csv`: BY-adjusted p-values from interaction model 
association testing. Row is phenotype, column is pairs of allelotypes shown to be
significantly associated with phenotype (if the cell is full). Only selected rows 
and columns that had significant associations (p<=0.05).

* `interaction_results_adjusted.csv`: LogORs, SEs, and p-values for the 3 terms
in the interaction GLM (term1, term2, interaction term) for each pairing of 
phenotype and  allelotype-pair.
