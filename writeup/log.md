### 20181205
* hwe test
* gender check revealed two problems with the sample table

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
miss 201939090006_R06C02 2018_CHOP_BAL_VEO_022
miss 201939090044_R03C01 2015_CHOP_MIC_BAL_FA071_SUB
miss 201939090044_R04C01 2015_CHOP_MIC_BAL_FAM106_SUB
miss 201939090044_R04C02 2015_CHOP_MIC_BAL_FAM018_SUB
miss 201939090044_R05C01 2015_CHOP_MIC_BAL_FAM035_SUB
miss 201939090044_R10C02 2015_CHOP_MIC_BAL_FAM092_sub
miss 201939090044_R11C02 2015_CHOP_MIC_BAL_FAM098_sub
miss 201939090044_R12C02 2015_CHOP__MIC_BAL_FAM061_SUB
miss 201939090076_R03C02 2015_CHOP_MIC_BAL_FAM042_SUB
miss 201939090076_R07C01 2015_CHOP_MIC_BAL_FAM)41_SUB
miss 201939090076_R08C02 2015_CHOP_MIC_BAL-FAM085_SUB
miss 201939090076_R09C01 2015_CHOP_MIC_BAL_FAM070_SUB
miss 201939090076_R10C01 2015_CHOP_MIC_BAL_FAM097_sub
miss 201939090114_R05C01 2915_CHOP_MIC_BAL_FAM176_SUB
miss 201939090114_R07C02 2015_CHOP_MIC_BAL_FAM109_SUB
miss 201939090114_R12C01 2015_CHOP_MIC_BAL_FAM164_SUB
miss 201939090156_R02C02 2015_CHOP_MIC_BAL_FAM87_SUB
miss 201939090156_R06C02 2015_CHOP_MIC_BAL_FAM166_SUB
miss 201939090166_R03C01 2015_CHOP_MIC_BAL_FAM94_SUB
miss 201939090193_R04C02 2015_CHOP_MIC_BAL_FAM076_SUB
miss 201939090193_R05C01 2015_CHOP_MIC_BAL_FAM101_SUB
miss 201939090193_R10C01 2015_CHOP_MIC_BAL_FAM0104_SUB
miss 201978470008_R02C02 2017_CHOP_BAL_VEO_007
miss 201978470008_R05C02 2017_CHOP_BAL_VEO_009
miss 201978470008_R07C01 2017_CHOP_BAL_VEO_010

### 20181120
* 173 samples after mds cut. 39 are cases and 134 are controls
* mds picks up 2 outlier case/control pairs that are very close
* 201939090050_R10C02 (missing 0.0009588) and 201978470008_R03C02 (missing 0.0008311) co cluster; drop 201939090050_R10C02
* 201939090076_R03C01 (missing 0.00138) and 201939090006_R11C02 (missing 0.0006981) co cluster; drop 201939090076_R03C01
* these are clusters of the same ppl. Remove the one w/ more missing targets

#### mtg
* keep XY and check sex. added chrX for hapmap and study
* Send missing samples, and % missing. Are they case/control?
* for indep, use indep-pairwise 0.2 updated.
* add ethnicity to mds plot w/ hapmap. done
* check pruned target #
* check hapmap/study overlap #
* rm two duplicate samples. done
* Send file z scores to check second cousins
* how many targets below 1%, 1-5%, and above 5%? what frac of total targets?
* hwe test. split case vs control/healthy. case is early and late ibd
* missing rates for case vs control
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

### 2018114
* setup plink pipeline

### 2018113
* gwas meeting
* what are we testing?
