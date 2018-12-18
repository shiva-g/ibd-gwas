### 20181218

#### Agenda
* [maf counts](tables/maf.md)
* prsice clumping defaults
* [prsice results](methods.md#polygenic-risk-score)
* [plink association results](methods.md#associations)
* [related samples](log.md#20181217) after mds w/ fixed indep snps

#### Notes

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
