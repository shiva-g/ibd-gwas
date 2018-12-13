### SNP filters
* Start with 683369 markers across 118 cases and 99 controls. 
* Use plink to remove 97585 monomorphic snps and markers not called for more than 5% of samples (geno 0.05). 
* Remove 1363 duplicated variants and 4560 positions that appeared more than once with plink.
* Remove 3 samples with more than 5% of markers not called (plink mind 0.05).
    *   IBD case 201939090179_R06C01          N    52791   557037  0.09477
    *   HC 201939090179_R04C01          N    62569   557037   0.1123
    *   IBD 201939090179_R05C01          N    66395   557037   0.1192
* For MDS, use independent markers (plink --indep-pairwise 50 5 0.2).

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

### QC
* sample missing table: /mnt/isilon/microbiome/perry/ibd-gwas/data/interim/missing_test/3groups.imiss
* hwe (--hardy) run after filtering samples with too many missing targets and removing chrX, but before indep. 
    * AFF: 1895 targets p<0.001. 7511 targets p<0.01 (might have real signal)
    * UNAFF controls: 2058 targets p<0.001. 7932 targets p<0.01
* gender check (--check-sex) run after indep filter shows two conflicts F of .08 and .06
* --het check: 0.08715, largest F inbreeding coefficient estimate; smallest -0.09699
* case vs control missing per target (--test-missing) run after filtering targets with high missing rates and removing chrX. 344 markers have p<0.05. 76 p<0.01. 7 p<0.001
* mafs after QC
    * Below 1%        50637
    * 1% to 5%        186738
    * Greater 5%      317792
    
### Imputation
* Michigan Imputation Server v1.0.4
* EUR sample subset after filtering missing samples (72 HC, 92 IBD)
* Reference Panel: hrc.r1.1.2016
* Phasing: eagle

### Polygenic risk score
* 227 of 232 base variants included (4 discarded b/c of mismatch). [Base snps](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4881818/). [Nature tables](https://www.nature.com/articles/ng.3359#supplementary-information)
* 164 people (98 male(s), 66 female(s)) 72 HC. 92 IBD
* 223 variants after clumping
* 1 region(s) with p-value less than 1e-5
* these results are inflated due to the overfitting
* [Barplot](plots/eur_BARPLOT_2018-12-11.png)
* [Hires plot](plots/eur_HIGH-RES_PLOT_2018-12-11.png)

#### PRSicse IBD vs HC
| Set |    Threshold   |    R2   |   P   |    Coefficient  |   Standard.Error | Num_SNP |
| ---- | ----- | ----- | ----- | ----- | ----- | ----- |
| Base  |  0.0001 | 0.195872  |      5.32139e-06  |   299.55 | 65.8106 | 222 |
| Base  |  0.0024 | 0.196663   |     5.09771e-06  |   300.708 | 65.9343 | 223 |

#### Summary IBD vs HC
| Phenotype  |     Set   |  Threshold   |    PRS.R2  | Full.R2 | Null.R2 | Prevalence  |    Coefficient |    Standard.Error | P    |   Num_SNP |
| -| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | --- |
| -    | Base  |  0.0024 | 0.196663   |     0.196663    |    0   |    -  |     300.708 | 65.9343 | 5.09771e-06  |   223 |

#### PRSicse IBD late vs early
