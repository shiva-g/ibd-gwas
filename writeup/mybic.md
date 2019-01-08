{% extends "base.html" %} {% load markdown_deux_tags %} {% block content %} {% markdown %}

### ibd_gwas
* [repo](https://github.com/samesense/ibd-gwas)
* [methods](https://github.com/samesense/ibd-gwas/blob/master/writeup/methods.md)

#### EUR association SNP tables
* [EUR IBD vs HC]({{SLINK}}/writeup/tables/all.eur.assoc.csv)
* [EUR early IBD vs HC]({{SLINK}}/writeup/tables/early.eur.assoc.csv)
* [EUR late IBD vs HC]({{SLINK}}/writeup/tables/late.eur.assoc.csv)
* [EUR late vs early IBD]({{SLINK}}/writeup/tables/ibd.eur.assoc.csv)

{% endmarkdown %} {% endblock %}
