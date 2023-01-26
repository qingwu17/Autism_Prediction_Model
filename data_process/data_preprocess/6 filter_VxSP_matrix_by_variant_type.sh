#!/bin/bash
#SBATCH -J filter_VxSP_matrix_by_variant_type
#SBATCH -n 8
#SBATCH -t 96:00:00
#SBATCH --mem=128G

# remove variants by index (including exonic and splicing variants, exclude synonymous and the rest)

# include exonic variants except synonymous variants
# add "d" to the end of each character
sed 's/$/d/' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.sed
sed -i "s/;$//g" /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.sed
# filter 
sed -f /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.sed \
/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix_cleaned.txt

# PTVs_MisA
# add "d" to the end of each character
sed 's/$/d/' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisA.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisA.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisA.sed
sed -i "s/;$//g" /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisA.sed
# filter 
sed -f /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisA.sed \
/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix_cleaned.PTVs_MisA.txt


# PTVs_MisAB
# add "d" to the end of each character
sed 's/$/d/' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisAB.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisAB.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisAB.sed
sed -i "s/;$//g" /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisAB.sed
# filter 
sed -f /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.PTVs_MisAB.sed \
/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix_cleaned.PTVs_MisAB.txt

# MisA
# add "d" to the end of each character
sed 's/$/d/' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisA.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisA.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisA.sed
sed -i "s/;$//g" /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisA.sed
# filter 
sed -f /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisA.sed \
/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix_cleaned.MisA.txt

# MisB
# add "d" to the end of each character
sed 's/$/d/' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisB.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisB.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisB.sed
sed -i "s/;$//g" /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisB.sed
# filter 
sed -f /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisB.sed \
/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix_cleaned.MisB.txt

# MisAB
# add "d" to the end of each character
sed 's/$/d/' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisAB.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisAB.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisAB.sed
sed -i "s/;$//g" /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisAB.sed
# filter 
sed -f /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.MisAB.sed \
/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix_cleaned.MisAB.txt

# MisC
# add "d" to the end of each character
sed 's/$/d/' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.non_impactful.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.non_impactful.sed
# substitute /n by ;
sed -z -i 's/\n/;/g' /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.non_impactful.sed
sed -i "s/;$//g" /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.non_impactful.sed
# filter 
sed -f /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants.idx_var2exclude.non_impactful.sed \
/path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix.txt \
> /path_to/wes1_wes2_combined.deepvariant.rare1pct_variants_by_sample_matrix_cleaned.non_impactful.txt

