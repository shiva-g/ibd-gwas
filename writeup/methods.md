### GSA+ SNP filters
* Start with 683369 markers across 118 cases and 99 controls. 
* Use plink to remove 97585 monomorphic snps and markers not called for more than 5% of samples (geno 0.05). 
* Remove 1363 duplicated variants and 4560 positions that appeared more than once with plink.
* Remove 3 samples with more than 5% of markers not called (plink mind 0.05).
    *   IBD case 201939090179_R06C01          N    52791   557037  0.09477
    *   HC 201939090179_R04C01          N    62569   557037   0.1123
    *   IBD 201939090179_R05C01          N    66395   557037   0.1192
* For MDS, use independent markers (plink --indep-pairwise 50 5 0.2).

### GSA SNP filters
* Start with 633220 markers across 266 cases and 719 controls.
* 545633 markers after variant filtering
* Remove 4 samples with more than 5% of markers not called (plink mind 0.05).
    * 1 VEO (15% of targets missing), 3 control

### HapMap samples
* Download Phase 2 HapMap as a PLINK fileset (CEU, YRI, and JPT_CHB filtered founders) from http://zzz.bwh.harvard.edu/plink/res.shtml.
* Use liftover to convert hg18 positions to hg19, discarding unmatched positions. Use plink to update genomic coordinates based on liftover results.
* Restrict to autosomes and chrX, and apply SNP filters below, ending with 209 people and 1553222 markers.

### MDS
* Combine HapMap and study plink files after filtering (see SNP filters). 
* Remove 10 variants with multiple positions, and 19 variants with 3+ alleles.
* Remove 21926 variants where the position matches, but allele1 or allele2 do not across the study and hapmap.
* Use plink to combine HapMap and study plink files with merge-mode 1 (concensus) and merge-equal-pos. Apply SNP filters.
* Use plink to make MDS coordinates from 158106 variants across 116 cases, 98 controls, and 209 HapMap samples.
* [IBD, HC, and hapmap plot](plots/hapmap_mds.png)
* [EUR IBC and HC plot](plots/eur_mds.png)

### GSA+ QC
* test fraction of filtered snps missed by sample
    * `sed -e 's/\s\+/\t/g' 3groups.imiss | cut -f 7 | sort -gr | head`
    * /mnt/isilon/microbiome/perry/ibd-gwas/data/interim/missing_test/3groups.imiss
    * 3 samples removed b/c more than 5% of filtered targets were missing
    * 12%, 11%, 9% of filtered targets were missing
* hwe (--hardy) run after filtering samples with too many missing targets and removing chrX, but before indep. 
    * 555083 variants
    * AFF: 1898 targets p<0.001. 7516 targets p<0.01 (might have real signal)
    * UNAFF controls: 1864 targets p<0.001. 7561 targets p<0.01
* gender check (--check-sex) run after indep filter shows two conflicts F of .08 and .06
* --het check: 0.08715, largest F inbreeding coefficient estimate; smallest -0.09699
* case vs control missing per target (--test-missing) run after filtering targets with high missing rates and removing chrX. 859 markers have p<0.05. 133 p<0.01. 7 p<0.001
* [plink --genome table](https://docs.google.com/spreadsheets/d/1CFsaf5nz1TcppBgqOd4VkWKGRO2t1xmxFcfQC6oYsjw/edit#gid=1911340057)
* [maf table](tables/maf.md)

### GSA QC
* test fraction of filtered snps missed by sample
    * /mnt/isilon/microbiome/perry/ibd-gwas/data/interim/missing_test/gsa.imiss
    * 4 samples removed b/c more than 5% of filtered targets were missing
    * 39%, 22%, 15%, 6% of filtered targets were missing
* hwe (--hardy) run after filtering samples with too many missing targets and removing chrX, but before indep.
    * 532312 variants
    * AFF: 4588 targets p<0.001. 13982 targets p<0.01 (might have real signal)
    * UNAFF controls: nan
* ~gender check (--check-sex) run after indep filter shows two conflicts F of .08 and .06~
* --het check: 0.365, largest F inbreeding coefficient estimate; smallest -0.15
* case vs control missing per target (--test-missing) run after filtering targets with high missing rates and removing chrX. 6981 markers have p<0.05. 3200 p<0.01. 1197 p<0.001
* [plink --genome table](https://docs.google.com/spreadsheets/d/1QK4bAMm4bZqctnldZjs5Jwbs1MiOqzWwflRetY1-RW0/edit#gid=301566719)
* ~[maf table](tables/maf.md)~
    
### Imputation
* Michigan Imputation Server v1.0.4
* EUR sample subset after filtering missing samples (64 HC, 91 IBD)
* Reference Panel: hrc.r1.1.2016
* Phasing: eagle
* Filter imputed results for imputed score R2>0.3 and maf>1%. Final snp count: 7,554,996

### Associations
* 64 HC. 91 IBD (53 late, 38 VEO)
* 7,554,996 imputed snps (imputed R2>0.3 and maf>1%)
* Use plink to find associations and qqman for plots
    * [HC vs IBD](plots/manhattan.all.png), [qq](plots/qq.all.png)
    * [HC vs VEO](plots/manhattan.early.png), [qq](plots/qq.early.png)
    * [HC vs late IBD](plots/manhattan.late.png), [qq](plots/qq.late.png)
    * [VEO vs late IBD](plots/manhattan.ibd_all.png), [qq](plots/qq.ibd_all.png)
* Use snptest for associations (-frequentist 1 -method score -hwe)    
* [See mybic `EUR association SNP tables` for SNP tables](http://mybic.chop.edu/labs/devoto_lab/ibd-gwas/)

### Polygenic risk score
* pvals from http://www.biosoft.hacettepe.edu.tr/easyROC/)
* 231 of 232 base variants. [Base snps](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4881818/). [Nature tables](https://www.nature.com/articles/ng.3359#supplementary-information)
    * 192 on array, 39 imputed
    * 16:50763781 indel not imputed or on array. The position is correct in the supplement, but the imputed have nothing within 30bp of it.

* 64 HC. 91 IBD (53 late, 38 VEO)
* [pval table](tables/prs.md)
* [IBD vs HC ROC](plots/all.eur.prs.roc.png) [IBD vs HC density](plots/all.eur.prs.density.png) [IBD vs HC quartiles](plots/all.eur.prs.quartiles.png) ROC pval: 1.596776e-09
* [VEO vs HC ROC](plots/early.eur.prs.roc.png) [VEO vs HC density](plots/early.eur.prs.density.png) [VEO vs HC quartiles](plots/early.eur.prs.quartiles.png) ROC pval: 0.000238243
* [Late IBD vs HC ROC](plots/late.eur.prs.roc.png) [Late IBD vs HC density](plots/late.eur.prs.density.png) [Late IBD vs HC quartiles](plots/late.eur.prs.quartiles.png) ROC pval: 6.77093e-09
* [Late IBD vs VEO ROC](plots/ibd_all.eur.prs.roc.png) [Late IBD vs VEO density](plots/ibd_all.eur.prs.density.png) [Late IBD vs VEO ROC quartiles](plots/ibd_all.eur.prs.quartiles.png) ROC pval: 	0.7775863
