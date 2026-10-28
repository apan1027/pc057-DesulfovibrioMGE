---
author:
- name: Jinlong Ru
- name: Cunli Pan
params:
  name: 170-figureS2-terl
subtitle: Redraw figS2_terl from frozen records
title: 170-figureS2-terl
---

This step reads the registered frozen inputs. It does not fit models,
run bootstrap resampling, call annotations or update the original
database. See the [output map](../metadata/result-code-map.tsv) and
[reproducibility scope](../docs/reproducibility.md).

## Tasks

``` r
pc057_draw_figure(step, "figS2_terl")
```

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                step status         output
    ## 1 170-figureS2-terl   pass figS2_terl.png

<figure>
<img src="data/170-figureS2-terl/figS2_terl.png"
alt="plot of chunk figure" />
<figcaption aria-hidden="true">plot of chunk figure</figcaption>
</figure>
