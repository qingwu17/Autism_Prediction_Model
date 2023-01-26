#!/bin/bash
#SBATCH -J gctb
#SBATCH -n 16
#SBATCH -t 36:00:00
#SBATCH --mem=64G

gctb_2.0_Linux/gctb \
--sbayes R \
--gwas-summary /path_to/daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.ma \
--mldm /path_to/gctb_2.0_Linux/gctb_2.0_tutorial/ldm/band_ukb_10k_hm3/ukb10k.mldm \
--pi 0.95,0.02,0.02,0.01 \
--gamma 0.0,0.01,0.1,1 \
--chain-length 10000 \
--burn-in 2000 \
--exclude-mhc \
--out-freq 10 \
--out /path_to/iPSYCH_PGC_comm_var/daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.out \
> /path_to/iPSYCH_PGC_comm_var/daner_iPSYCH-DBS_Aut_PGC1_Eur_ph3.log 2>&1
