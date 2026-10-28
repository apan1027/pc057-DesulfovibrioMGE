---
author:
- name: Jinlong Ru
- name: Cunli Pan
params:
  name: 140-figure4-neighborhood
subtitle: Redraw fig4_neighborhood from frozen records
title: 140-figure4-neighborhood
---

This step reads the registered frozen inputs. It does not fit models,
run bootstrap resampling, call annotations or update the original
database. See the [output map](../metadata/result-code-map.tsv) and
[reproducibility scope](../docs/reproducibility.md).

## Tasks

``` r
pc057_draw_figure(step, "fig4_neighborhood")
```

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                       step status                output
    ## 1 140-figure4-neighborhood   pass fig4_neighborhood.pdf
    ## 2 140-figure4-neighborhood   pass fig4_neighborhood.png

<figure>
<img src="data/140-figure4-neighborhood/fig4_neighborhood.png"
alt="plot of chunk figure" />
<figcaption aria-hidden="true">plot of chunk figure</figcaption>
</figure>
