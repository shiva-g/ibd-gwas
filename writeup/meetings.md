### 20190108

#### Agenda
* [mybic docs](http://mybic.chop.edu/labs/devoto_lab/ibd-gwas/)
* [EUR plink & snptest annotated assoc results](http://mybic.chop.edu/labs/devoto_lab/ibd-gwas/)
* [QQ plots for plink assoc test](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md#associations)
* [PRSice quartile plots](https://github.com/samesense/ibd-gwas/blob/master/writeup/tables/prs.md)
* [AUC for PRSice ROC curves](https://github.com/samesense/ibd-gwas/blob/master/writeup/tables/prs.md)
* Michigan alternative for local compute imputation?

#### Working
* impute everyone w/ 1kg ref and run prsice w/ multi ethinic panel base snps
* pval for rocs

#### Notes

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
* ~~auc~~ and pval for rocs
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
