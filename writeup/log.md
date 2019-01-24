### 20190123
* Missing race for one gsa sample
* Fixed race sample name issues

### 20190122
* gsa mds done, but sample info wrong

### 20190121
* after gsa, hapmap merge, no snp has less than 5% of samples missing
* gsa bim files are 0/1! Not nucs. Fixed
* Next fix race in /mnt/isilon/microbiome/perry/ibd-gwas/data/interim/mds_dat/ibd_gsa_hapmap.dat

### 20190116
* merging new gsa data is not simple b/c the snp IDs disagree
* need to adjust mds rules to account for different dataset

### 20190115
* NOD2 variant [rs2066844](https://www.snpedia.com/index.php/Rs2066844) at 16:50745926. T/A is the risk allele.
* CARD9 9:139266405 [rs10781499](https://www.ebi.ac.uk/gwas/variants/rs10781499) T/A is the risk allele.

### 20190114
* IBD supplement does not have multi ethnic pvals, so I'll use eur for now
* Using PRS w/ 2 PCs looks very similar to PRS w/o PCs for total populations and EUR base SNPs

### 20190109
* michigan server is up. Running all pop in queue position 4. dl results
* FUT2 chr19:49206674 [rs601338](https://www.ncbi.nlm.nih.gov/snp/rs601338) G>A ref is G. What does each allele do? [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4992076/) This SNP is not on the chip, but [rs516246](https://www.ncbi.nlm.nih.gov/snp/rs516246) is (chr19:49206172). [Individuals with the homozygote A/A genotype are defined as non-secretors. Non-secretors, who are homozygous for the loss-of-function alleles of FUT2 gene (sese), have increased susceptibility to Crohn's disease](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4992076/) T or A is the risk allele. I can look up cases and controls maf from the snptest results. How do I split the cohort based on a SNP? For snptest: T is 0.48897299999999994 cases, 0.492188 HC. I can use plink extract to get the samples w/ this SNP. Plink assoc and snptest have flipped A1/A2.

### 20190108
* meeting
* not enough docker space to grab 1kg for imputation

### 20190107
* summarize snptest eur
* local impute tpop
    * `docker exec -t -i mis-docker cloudgene run imputationserver --files ${TEST_DATA} --refpanel apps@hapmap2 --conf /etc/hadoop/conf`
    * tpop PCA
    * tpop snptest
* mybic site - ask shiping for help
    * add snp results

### 20190103
* run snptest conditioning on covariates for pop structure using [cov_names cov1 and cov2](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#conditional_tests) when all sample imputation is done

### 20190102
* rm? reported white, but appears black: HC (201939090044_R12C01, 2015_CHOP_MIC_BAL_FAM029_SUB). No keep this.
* mk tpop vcfs for imputation. imputation server down

### 20181221
* what genome version for nature tables? hg19
* [nature tables messed up](https://mail.google.com/mail/u/0/#sent/QgrcJHsBmtxPmDHJwQbCGTJqzHkgZxlBDCg)

### 20181218
* process latest imputed results
* mafs eur, after impute, all samples

### 20181217
* new mds cut
    * odd samples: 201939090193_R10C02, 201939090179_R07C02
    
| FID1 | IID1 | FID2 | IID2 | RT | EZ | Z0 | Z1 | Z2 | PI_HAT | PHE | DST | PPC | RATIO |
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
| 0 | 201939090179_R07C02 | 0 | 201939090193_R10C02 | OT | 0 | 0.7159 | 0.2127 | 0.0714 | 0.1777 | -1 | 0.772610 | 1.0000 | 2.7478 |

* missing: 0.0009115 201939090193_R10C02, 0.0007926 201939090179_R07C02. Rm 201939090193_R10C02
* Next clustering produces two paired groups:
* 201939090050_R02C01 (rm), 201939090050_R01C02
    
| FID1 | IID1 | FID2 | IID2 | RT | EZ | Z0 | Z1 | Z2 | PI_HAT | PHE | DST | PPC | RATIO |
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
| 0 | 201939090050_R01C02 | 0 | 201939090050_R02C01 | OT | 0 | 0.7032 | 0.2377 | 0.0591 | 0.1780 | -1 | 0.771764 | 1.0000 | 2.5876 |

* 201939090050_R04C02 (rm), 201939090114_R10C01

| FID1 | IID1 | FID2 | IID2 | RT | EZ | Z0 | Z1 | Z2 | PI_HAT | PHE | DST | PPC | RATIO |
|--|--|--|--|--|--|--|--|--|--|--|--|--|--|
| 0 | 201939090050_R04C02 | 0 | 201939090114_R10C01 | OT | 0 | 0.7553 | 0.1461 |  0.0986 | 0.1717 | -1 | 0.773065 | 0.9999 | 2.5153 |
    
* manhatn plt: why nans in assoc file?
* annotate gene/pathway

### 20181214
* does plink bfile->vcf drop vars? no
    * 157166 chr22 r2 limit vcf
    * 157166 chr22 bim

### 20181213
* try impute server: https://imputationserver.readthedocs.io/en/latest/docker/
* `docker run --user $(id -u) -d -p 8080:80 -e DOCKER_CORES="4" -v /home/evansj/server/:/data/ --name mis-docker genepi/imputationserver`
* it runs on refosco
* local server does not allow HC ref data :(
* plink association tests
* prs for other groups
* ibd enhancers here? https://www.nature.com/articles/s41467-018-03766-z#additional-information

### 20181212
* MAF calcs before indep, on QC snps
    * code changed
    * indep filter does nothing
    * I was not applying the snp list produced by indep--pairwise
    * hapmap clustering w/ real indep snps kicks out 5 eur samples
* empirical prsice pval
* het check
* R2 in vcf limit of 0.3 and 1% maf on imputation results
* sample mixups
    * HC 109 is strange: In the sample table from 11/28, the sample is in the Withdrawn tab b/c there are no samples. In the 12/10 table, the sample is moved to the Has GWAS tab. I'll drop this sample since it is a control.
    * IBD case 087 is reported as white/black and clusters with the african group. It isn't included for the current EUR analysis, so the gender discrepancy won't matter until we do an AFR analysis.

### 20181211
* Monomorphic sites: 196 on chr22. I thought I removed these. I bet these are monomorphic for just eur.
* new impute run done
* mtg
    * ~~check het.~~ Do flipped gender samples have bad het scores? AFR should have more het than EUR
    * ~~MAF calcs before indep, on QC snps~~
    * ~~snp count before imputing~~, and # of imputed snps
    * new mafs after imputation
    * ~~R2 in vcf limit of 0.3 and 1% maf on imputation results~~
    * ~~what is clumping r2 for prsice?~~ what are the mismatched snps?
        * --clump-r2 0.1
        * --clump-kb 250 (therefore a 500kb window with the index SNP at the center)
    * ~~empirical prsice pval~~
    * ~~rm prsice plots~~
    * ~~use prsice scores in best to look at hc, ibd, late, early dists and roc~~
    * ~~veo vs older prs~~
    * ~~plink assoc veo vs old, veo vs hc, ibd vs hc, all vs hc~~
 
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
