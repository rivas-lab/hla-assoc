# hla-assoc

## Rounding for dosage in .ped file

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
