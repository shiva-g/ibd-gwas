### 20190215

#### Updates
* Removed rs2226628 from PRS
* Compare ped/adult risk alleles (if alleles match, does odds ratio match?)
    * problem indel at adult snp 16:50763781 A1=D, A2=I (present before imputation, but not after)
    * EUR: 141 agree, 86 disagree, 4 unresolved
    * All pops: 147 match, 76 disagree, 4 unresolved
    * plink vs snptest:
        * https://docs.google.com/spreadsheets/d/1LFdhrfXSHDtJ3Ed3tt7ynj9LL4LY3LMvuiCwaEBFpYE/edit#gid=551803840
* Pathway PRS (on going)
    * Building pathways from [Khor 2011](https://www.ncbi.nlm.nih.gov/pubmed/21677747) gene lists and associated GO terms
    * 64/231 SNPs matched to a pathway
    * 135/231 associated w/ GO term
    * 64/135 possible pathway assignments made
    * [EUR](https://github.com/samesense/ibd-gwas/blob/master/writeup/tables/prs.eur.md) sorted by p-value
    * [All pops](https://github.com/samesense/ibd-gwas/blob/master/writeup/tables/prs.tpop.md) sorted by p-value
    * [SNP pathway assignments](http://mybic.chop.edu/slink/devoto_lab/ibd-gwas/writeup/tables/adult.all.tpop.assoc.csv)
* Association tests for adult genes (ongoing)
* nealelab data on dropbox

#### Notes
* ~bug in pathway snps - I joing on plink A1/A2 instead of adult A1/A2~
* 5 ped specific SNPs (not seen in adult) from Hakon paper
* merge pathways to make larger groups - get gene sets from Noor
* Check [gene ls](https://mail.google.com/mail/u/0/#inbox/FMfcgxwBVgqHjxqmbzzPCBJFZSZtPpdC) against nealelab results
* don't restrict SNP OR adult/ped agreement for PRS for you want to test the utility of the adult SNPs as a set

### 20190201

#### Overview
* GSA+: 115 are cases and 94 are controls
* GSA: 265 are cases (41 on GSA+) and 702 are controls (unknown control overlap)
    * 24A female vs 0 (F=0.5) female in both gender tables
    * 1070A female vs male (F=1) female in both gender tables; male by genotype
* 231 adult IBD SNPs used for polygenic risk score to rank GSA+. [AUC 0.73982 all populations](https://github.com/samesense/ibd-gwas/blob/master/writeup/plots/all.tpop.prs.roc.png). [All population PRS density](https://github.com/samesense/ibd-gwas/blob/master/writeup/plots/all.tpop.prs.density.png)
* GSA+ SNP associations
    * [231 adult ranked by snptest pval for all populations](https://github.com/samesense/ibd-gwas/blob/master/writeup/tables/adult.all.tpop.assoc.csv)
    * [All SNPs below p 1e-5; all populations](https://github.com/samesense/ibd-gwas/blob/master/writeup/tables/ped.all.tpop.assoc.csv)
    * [All SNPs plot EUR](https://github.com/samesense/ibd-gwas/blob/master/writeup/plots/manhattan.all.png)
* Split GSA+ samples by PRS quantile, FUT2, CARD9, NOD2 risk allele

### 20190129

#### Updates
* [VEO GSA processing](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md#veo-gsa-snp-filters)
* [VEO GSA QC](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md#veo-gsa-qc)
    * [GSA genome table](https://docs.google.com/spreadsheets/d/1QK4bAMm4bZqctnldZjs5Jwbs1MiOqzWwflRetY1-RW0/edit#gid=301566719)
        * 2 VEO identical to contols
        * 2 related VEO
        * 1 VEO related to control
        * Related controls
    * [GSA+ genome table](https://docs.google.com/spreadsheets/d/1CFsaf5nz1TcppBgqOd4VkWKGRO2t1xmxFcfQC6oYsjw/edit#gid=1911340057)
    * .124 PI_HAT cutoff?
* [GSA MDS](https://github.com/samesense/ibd-gwas/blob/master/writeup/plots/hapmap_mds.gsa.png)  
* IBD GSA+ subject ID mixup
    * Subject 261 labels two different people
    * Subject 215 labels two different people
* Meetings?
    * Tues 3:30 Y
    * Dynamic VEO microbiome 10:15?
    * IBD Tues 11?
    * IBD Fri 11am Y
    * Large dynamic group?
* Judith IBD pathways; split samples by pathway?
* [gsa vs gsa+ samples](https://docs.google.com/spreadsheets/d/1T5TmfiabuY-EuKc16Dmynmtao7TAAqAqU-z_hBAHQnI/edit#gid=186652860)
    * 41 VEO duplicates
    * 0 control duplicates by ID

#### Working
* Resolving sample issues
* Apply filters discussed during meeting
* Clean microbiome table for Yue

#### Notes
* Contact Ying about control samples found as cases
* Ask about 24A cannot tell gender next mtg
* Ask about 1070A female w/ male genotype (found in last analysis)

### 20190115

#### Updates
* Multi-ethnic analysis (94 HC, 116 IBD) (43 VEO, 73 late IBD)
    * PRS: Binary calls, not bgen, snptest: bgen
    * First 2 PCs
    * Use adult EUR 231 base snps
    * [Multi-ethnic polygenic risk score table](tables/prs.tpop.md)
    * [Multi-ethnic polygenic risk score quantiles](plots/all.tpop.prs.quantiles.png)
    * [Multi-ethnic polygenic risk score roc](plots/all.tpop.prs.roc.png)
    * [Multi-ethnic snptest for adult SNPs](tables/adult.all.tpop.assoc.csv) ([EUR](tables/adult.all.eur.assoc.csv))
    * [Multi-ethnic snptest](tables/ped.all.tpop.assoc.csv) ([EUR](tables/ped.all.eur.assoc.csv))
* [EUR polygenic risk score quantiles](plots/all.eur.prs.quantiles.png)
* Have UK biobank data access

#### Working
* ~for fut2, card9, and nod2, how do adult risk allele vs non-risk split up patients?~
* assign genes to pathways; how do patients split by pathways using risk alleles?
* do adult and ped effect agree?

#### Notes
* [more data](https://mail.google.com/mail/u/0/#inbox/FMfcgxwBVDKqJbCSJrRCzNkLHvWCjNGc) to run by itself. Check overlap with existing samples.

### 20190108

#### Agenda
* [mybic docs](http://mybic.chop.edu/labs/devoto_lab/ibd-gwas/)
* [EUR plink annotated assoc results](http://mybic.chop.edu/labs/devoto_lab/ibd-gwas/). [ex table all ibd vs hc](tables/all.eur.assoc.csv)
* [QQ plots for plink assoc test](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md#associations)
* [PRSice quartile plots](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md#polygenic-risk-score)
* [AUC for PRSice ROC curves](https://github.com/samesense/ibd-gwas/blob/master/writeup/tables/prs.md)
* [pvals for rocs](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md#polygenic-risk-score)
* Michigan alternative for local compute imputation?
* [Welcome trust data form](https://www.dropbox.com/s/u60f4i2uhh7jtmc/Screenshot%202019-01-08%2014.15.14.png?dl=0)
* UK biobank registration under review

#### Working
* ~impute everyone w/ 1kg ref and run prsice w/ multi ethinic panel base snps~
* Welcome trust 1st edition ibd adult and Uk biobank ibd adult data
* automate roc pvals
* ~add snptest results to eur plink assoc table~

#### Notes
* ~Fix agenda links~
* ~Fix devoto mybic group~
* ~Add judith to mybic~
* ~Add ibd group colors to quantile plots~
* ~Send ibd quantile groups to yue with subject IDs~
* Report of analysis this far
* use easyROC cuts for splitting IBD vs HC to split IBD samples
* utilize impute2
* ~for fut2, card9, and nod2, how do how do adult risk allele vs non-risk split up patients?~
* assign genes to pathways; how do patients split by pathways using risk alleles?
* ~split snp assoc results and adult alleles into two tables~
* do adult and ped effect agree?
* label risk alleles and be able to split patients by hom, any risk, no risk

### 20181218

#### Agenda
* [maf counts](tables/maf.md)
* prsice clumping defaults
* [prsice results](methods.md#polygenic-risk-score)
* [plink association results](methods.md#associations)
* [related samples](log.md#20181217) after mds w/ fixed indep snps
* qqplot?

#### Working
* annotate plink association positions w/ ~genes~/pathways. [fuma?](https://www.nature.com/articles/s41467-017-01261-5)
* ~~nvestigate prsice mismatched snps~~

#### Notes
* ~impute everyone w/ 1kg ref and run prsice w/ multi ethinic panel base snps~
* ~plot/list prsice quartiles for splitting samples by quartile~
* ~what are prsice snps scores for assocation test?~ plot them
* ~~auc~~ and ~pval for rocs~
* ~why lost base snps in prsice?~
* ~why does prsice list different snps counts for each test?~
* ~~rm clumping step from prsice so all snps are kept~~
* ~use snptest w/ snp probs and 2 PCs and logistic regression instead of plink~
* ~~mk qqplot from association results~~

### 20181211

#### Notes
* ~~check het.~~ Do flipped gender samples have bad het scores? AFR should have more het than EUR
* ~~MAF calcs before indep, on QC snps~~
* ~~snp count before imputing~~, ~~and # of imputed snps~~
* ~~new mafs after imputation~~
* ~~R2 in vcf limit of 0.3 and 1% maf on imputation results~~
* ~~what is clumping r2 for prsice?~~ what are the mismatched snps?
    * --clump-r2 0.1
    * --clump-kb 250 (therefore a 500kb window with the index SNP at the center)
* ~~empirical prsice pval~~
* ~~rm prsice plots~~
* ~~use prsice scores in best to look at hc, ibd, late, early dists and roc~~
* ~~veo vs older prs~~
* ~~plink assoc veo vs old, veo vs hc, ibd vs hc, all vs hc~~
