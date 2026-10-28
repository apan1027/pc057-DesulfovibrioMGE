---
author:
- name: Jinlong Ru
- name: Cunli Pan
params:
  name: 150-figure5-outgroup
subtitle: Redraw fig5_outgroup from frozen records
title: 150-figure5-outgroup
---

This step reads the registered frozen inputs. It does not fit models,
run bootstrap resampling, call annotations or update the original
database. See the [output map](../metadata/result-code-map.tsv) and
[reproducibility scope](../docs/reproducibility.md).

## Tasks

``` r
pc057_draw_figure(step, "fig5_outgroup")
```

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                   step status            output
    ## 1 150-figure5-outgroup   pass fig5_outgroup.pdf
    ## 2 150-figure5-outgroup   pass fig5_outgroup.png

<figure>
<img src="data/150-figure5-outgroup/fig5_outgroup.png"
alt="plot of chunk figure" />
<figcaption aria-hidden="true">plot of chunk figure</figcaption>
</figure>
