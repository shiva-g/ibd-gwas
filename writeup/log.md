### 20181120
* 173 samples after mds cut. 39 are cases and 134 are controls
* mds picks up 2 outlier case/control pairs that are very close
* 201939090050_R10C02 and 201978470008_R03C02 co cluster
* 201939090076_R03C01 and 201939090006_R11C02 co cluster
* these are clusters of the same ppl. Remove the one w/ more missing targets

#### mtg
* keep XY and check sex
* Send missing samples, and % missing. Are they case/control?
* for indep, use indep-pairwise 0.2
* add ethnicity to mds plot w/ hapmap
* check pruned target #
* check hapmap/study overlap #
* rm two duplicate samples
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
