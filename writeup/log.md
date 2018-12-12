### 20181212
* MAF calcs before indep, on QC snps
    * code changed
    * indep filter does nothing

### 20181211
* Monomorphic sites: 196 on chr22. I thought I removed these. I bet these are monomorphic for just eur.
* new impute run done
* mtg
    * check -het. Do flipped gender samples have bad het scores? AFR should have more het than EUR
    * MAF calcs before indep, on QC snps
    * snp count before imputing, and # of imputed snps
    * new mafs after imputation
    * R2 in vcf limit of 0.3 and 1% maf on imputation results
    * what is clumping r2 for prsice? what are the mismatched snps?
    * empirical prsice pval
    * rm prsice plots
    * use prsice scores in best  to look at hc, ibd, late, early dists and roc
    * veo vs older prs
    * plink assoc veo vs old, veo vs hc, ibd vs hc, all vs hc
 
### 20181210
* rm another sample
* switch from 0/1 format to nuc
* two more gender issues
* in impute server queue all day

### 20181207
* have new sample data

#### code review notes
* quay - plink and prsice
* plink logs vs perry logs
* iterative plink
* snakemake wrappers
* imputation server

### 20181206
* testing prsice - overlap is done by snp id, and they do not match. use --keep-ambig to turn this off

### 20181205
* hwe test. do not use x
* gender check revealed two problems with the sample table
* missing check: do not use x

### 20181204
* case (all ibd) vs control missing test

### 20181203
* [new sample table](https://mail.google.com/mail/u/0/#inbox/FMfcgxvzLrJTTbPBRKbqKrlxdgCqjQKh)
* all ibd samples (early and late) are now cases
* 01939090166_R10C01 and 201939090050_R05C02 ONC
* missing sample info? 
    * 201939090166_R03C01: 2015_CHOP_MIC_BAL_FAM94_SUB
* These are marked as onc. Should I treat them as healthy?
    * onc?? 185 201939090050_R03C02 ONC 2015_CHOP_MIC_BAL_FAM126_SUB
    * onc?? 185 201939090050_R05C02 ONC 2015_CHOP_MIC_BAL_FAM181_SUB
    * onc?? 185 201939090166_R01C01 ONC 2015_CHOP_MIC_BAL_FAM200_SUB
    * onc?? 185 201939090166_R10C01 ONC 2015_CHOP_MIC_BAL_FAM170_SUB
* I'm missing data for these:
    * miss 201939090006_R06C02 2018_CHOP_BAL_VEO_022:discard
    * miss 201939090044_R03C01 2015_CHOP_MIC_BAL_FA071_SUB:2015_CHOP_MIC_BAL_FAM071_SUB
    * miss 201939090044_R04C01 2015_CHOP_MIC_BAL_FAM106_SUB:discard
    * miss 201939090044_R04C02 2015_CHOP_MIC_BAL_FAM018_SUB:discard
    * miss 201939090044_R05C01 2015_CHOP_MIC_BAL_FAM035_SUB:ok
    * miss 201939090044_R10C02 2015_CHOP_MIC_BAL_FAM092_sub:2015_CHOP_MIC_BAL_FAM092_SUB
    * miss 201939090044_R11C02 2015_CHOP_MIC_BAL_FAM098_sub:2015_CHOP_MIC_BAL_FAM098_SUB
    * miss 201939090044_R12C02 2015_CHOP__MIC_BAL_FAM061_SUB:2015_CHOP_MIC_BAL_FAM061_SUB
    * miss 201939090076_R03C02 2015_CHOP_MIC_BAL_FAM042_SUB:ok
    * miss 201939090076_R07C01 2015_CHOP_MIC_BAL_FAM)41_SUB:2015_CHOP_MIC_BAL_FAM041_SUB
    * miss 201939090076_R08C02 2015_CHOP_MIC_BAL-FAM085_SUB:2015_CHOP_MIC_BAL_FAM085_SUB
    * miss 201939090076_R09C01 2015_CHOP_MIC_BAL_FAM070_SUB:ok
    * miss 201939090076_R10C01 2015_CHOP_MIC_BAL_FAM097_sub:2015_CHOP_MIC_BAL_FAM097_SUB
    * miss 201939090114_R05C01 2915_CHOP_MIC_BAL_FAM176_SUB:2015_CHOP_MIC_BAL_FAM176_SUB
    * miss 201939090114_R07C02 2015_CHOP_MIC_BAL_FAM109_SUB:ok
    * miss 201939090114_R12C01 2015_CHOP_MIC_BAL_FAM164_SUB:ok
    * miss 201939090156_R02C02 2015_CHOP_MIC_BAL_FAM87_SUB:2015_CHOP_MIC_BAL_FAM087_SUB
    * miss 201939090156_R06C02 2015_CHOP_MIC_BAL_FAM166_SUB:ok
    * miss 201939090166_R03C01 2015_CHOP_MIC_BAL_FAM94_SUB:2015_CHOP_MIC_BAL_FAM094_SUB
    * miss 201939090193_R04C02 2015_CHOP_MIC_BAL_FAM076_SUB:discard
    * miss 201939090193_R05C01 2015_CHOP_MIC_BAL_FAM101_SUB:discard
    * miss 201939090193_R10C01 2015_CHOP_MIC_BAL_FAM0104_SUB:2015_CHOP_MIC_BAL_FAM104_SUB
    * miss 201978470008_R02C02 2017_CHOP_BAL_VEO_007:2018_CHOP_BAL_VEO_007
    * miss 201978470008_R05C02 2017_CHOP_BAL_VEO_009:2018_CHOP_BAL_VEO_009
    * miss 201978470008_R07C01 2017_CHOP_BAL_VEO_010:2018_CHOP_BAL_VEO_010

### 20181120
* 173 samples after mds cut. 39 are cases and 134 are controls
* mds picks up 2 outlier case/control pairs that are very close
* 201939090050_R10C02 (missing 0.0009588) and 201978470008_R03C02 (missing 0.0008311) co cluster; drop 201939090050_R10C02
* 201939090076_R03C01 (missing 0.00138) and 201939090006_R11C02 (missing 0.0006981) co cluster; drop 201939090076_R03C01
* these are clusters of the same ppl. Remove the one w/ more missing targets

#### mtg
* keep XY and check sex. added chrX for hapmap and study. done
* Send missing samples, and % missing. Are they case/control? done.
* for indep, use indep-pairwise 0.2 updated.
* add ethnicity to mds plot w/ hapmap. done
* check pruned target #
* check hapmap/study overlap #
* rm two duplicate samples. done
* Send file z scores to check second cousins
* how many targets below 1%, 1-5%, and above 5%? what frac of total targets? done
* hwe test. split case vs control/healthy. case is early and late ibd. done
* missing rates for case vs control. done
* imputation w/ michigain server. use X. use haplotype consortium
* polygenic risk score w/ snps from paper 3rd sheet: all loci eur
* next mtg dec 11

### 20181119
* methods
* use hwe for a filter? http://bio3.giga.ulg.ac.be/chaichoompu/userfiles/downloads/2016/GBIO0009/L5/Plink.pdf slide 19

### 20181116
* MDS shows most are CEU

### 20181115
* Use hapmap plink files for ancestry MDS http://zzz.bwh.harvard.edu/plink/res.shtml

### 20181114
* setup plink pipeline

### 20181113
* gwas meeting
* what are we testing?
