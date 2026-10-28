---
author:
- name: Jinlong Ru
- name: Cunli Pan
params:
  name: 130-figure3-enrichment
subtitle: Redraw fig3_enrichment from frozen records
title: 130-figure3-enrichment
---

This step reads the registered frozen inputs. It does not fit models,
run bootstrap resampling, call annotations or update the original
database. See the [output map](../metadata/result-code-map.tsv) and
[reproducibility scope](../docs/reproducibility.md).

## Tasks

``` r
pc057_draw_figure(step, "fig3_enrichment")
```

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                     step status              output
    ## 1 130-figure3-enrichment   pass fig3_enrichment.pdf
    ## 2 130-figure3-enrichment   pass fig3_enrichment.png

<figure>
<img src="data/130-figure3-enrichment/fig3_enrichment.png"
alt="plot of chunk figure" />
<figcaption aria-hidden="true">plot of chunk figure</figcaption>
</figure>
