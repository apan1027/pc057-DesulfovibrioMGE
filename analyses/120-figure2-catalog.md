---
author:
- name: Jinlong Ru
- name: Cunli Pan
params:
  name: 120-figure2-catalog
subtitle: Redraw fig2_catalog from frozen records
title: 120-figure2-catalog
---

This step reads the registered frozen inputs. It does not fit models,
run bootstrap resampling, call annotations or update the original
database. See the [output map](../metadata/result-code-map.tsv) and
[reproducibility scope](../docs/reproducibility.md).

## Tasks

``` r
pc057_draw_figure(step, "fig2_catalog")
```

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                  step status           output
    ## 1 120-figure2-catalog   pass fig2_catalog.pdf
    ## 2 120-figure2-catalog   pass fig2_catalog.png

<figure>
<img src="data/120-figure2-catalog/fig2_catalog.png"
alt="plot of chunk figure" />
<figcaption aria-hidden="true">plot of chunk figure</figcaption>
</figure>
