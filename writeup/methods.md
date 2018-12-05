### SNP filters
* Start with 683369 markers across 44 cases and 185 controls. Sex is not provided yet, so label all male. 
* Merge files with plink, and restrict to autosomes.
* Use plink to remove non-polymorphic markers (min-ac 1) and markers not called for more than 5% of samples (geno 0.05). 
* Remove 2728 duplicated variants and 4471 positions that appeared more than once with plink.
* Remove 3 samples with more than 5% of markers not called (plink mind 0.05).
* For MDS, use independent markers (plink indep 50 5 2).

### HapMap samples
* Download Phase 2 HapMap as a PLINK fileset (CEU, YRI, and JPT_CHB filtered founders) from http://zzz.bwh.harvard.edu/plink/res.shtml.
* Use liftover to convert hg18 positions to hg19, discarding unmatched positions. Use plink to update genomic coordinates based on liftover results.
* Restrict to autosomes, and apply SNP filters below, ending with 209 people and 1553222 markers.

### MDS
* Combine HapMap and study plink files after filtering (see SNP filters). 
* Remove 10 variants with multipole positions, and 19 variants with 3+ alleles.
* Remove 20668 variants where the position matches, but allele1 or allele2 do not across the study and hapmap.
* Use plink to combine HapMap and study plink files with merge-mode 1 (concensus) and merge-equal-pos. Apply SNP filters.
* Use plink to make MDS coordinates from 153443 variants across 44 cases, 182 controls, and 209 HapMap samples.

### QC
* hwe (--hardy) run after filtering samples with too many missing targets and removing chrX, but before indep. 
    * AFF: 1927 targets p<0.001. 7543 targets p<0.01
    * UNAFF: 2432 targets p<0.001. 8823 targets p<0.01
* gender check (--check-sex) run after indep filter shows no conflicts
* case vs control missing per target (--test-missing) run after filtering targets with high missing rates and removing chrX. 344 markers have p<0.05. 76 p<0.01. 7 p<0.001
