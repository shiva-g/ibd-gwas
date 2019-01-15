### 20190115

### Updates
* Multi-ethic analysis (94 HC, 116 IBD) (43 VEO, 73 late IBD)
    * Binary calls, not bgen
    * First 2 PCs
    * [Multi-ethnic polygenic risk score table](tables/prs.tpop.md)
    * [Multi-ethnic polygenic risk score quantiles](plots/all.tpop.prs.quantiles.png)
    * [Multi-ethnic polygenic risk score roc](plots/all.tpop.prs.roc.png)
    * [Multi-ethnic snptest for adult SNPs](tables/adult.all.tpop.assoc.csv) ([EUR](tables/adult.all.eur.assoc.csv))
    * [Multi-ethnic snptest](tables/ped.all.tpop.assoc.csv) ([EUR](tables/ped.all.eur.assoc.csv))
* [EUR polygenic risk score quantiles](plots/all.eur.prs.quantiles.png)

### Working
* for fut2, card9, and nod2, how do how do adult risk allele vs non-risk split up patients?
* assign genes to pathways; how do patients split by pathways using risk alleles?
* do adult and ped effect agree?

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
* impute everyone w/ 1kg ref and run prsice w/ multi ethinic panel base snps
* Welcome trust 1st edition ibd adult and Uk biobank ibd adult data
* automate roc pvals
* add snptest results to eur plink assoc table

#### Notes
* ~Fix agenda links~
* ~Fix devoto mybic group~
* ~Add judith to mybic~
* ~Add ibd group colors to quantile plots~
* ~Send ibd quantile groups to yue with subject IDs~
* Report of analysis this far
* use easyROC cuts for splitting IBD vs HC to split IBD samples
* utilize impute2
* for fut2, card9, and nod2, how do how do adult risk allele vs non-risk split up patients?
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
* impute everyone w/ 1kg ref and run prsice w/ multi ethinic panel base snps
* ~plot/list prsice quartiles for splitting samples by quartile~
* ~what are prsice snps scores for assocation test?~ plot them
* ~~auc~~ and ~pval for rocs~
* ~why lost base snps in prsice?~
* ~why does prsice list different snps counts for each test?~
* ~~rm clumping step from prsice so all snps are kept~~
* use snptest w/ snp probs and 2 PCs and logistic regression instead of plink
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
