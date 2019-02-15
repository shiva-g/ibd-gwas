{% extends "base.html" %} {% load markdown_deux_tags %} {% block content %} {% markdown %}

### ibd_gwas
* [repo](https://github.com/samesense/ibd-gwas)
* [methods](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md)

#### EUR association SNP tables
* SNPs were used for PRS, or are in adult genes
* [EUR IBD vs HC]({{SLINK}}/writeup/tables/ped.all.eur.assoc.csv)
* [EUR early IBD vs HC]({{SLINK}}/writeup/tables/ped.early.eur.assoc.csv)
* [EUR late IBD vs HC]({{SLINK}}/writeup/tables/ped.late.eur.assoc.csv)
* [EUR late vs early IBD]({{SLINK}}/writeup/tables/ped.ibd.eur.assoc.csv)

#### Total population associations for adult SNPs
* SNPs were used in PRS
* [Tot pop IBD vs HC]({{SLINK}}/writeup/tables/adult.all.tpop.assoc.csv)
* [Tot pop early IBD vs HC]({{SLINK}}/writeup/tables/adult.early.eur.assoc.csv)
* [Tot pop late IBD vs HC]({{SLINK}}/writeup/tables/adult.late.eur.assoc.csv)
*
{% endmarkdown %} {% endblock %}
