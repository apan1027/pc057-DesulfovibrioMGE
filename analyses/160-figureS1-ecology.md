---
author:
- name: Jinlong Ru
- name: Cunli Pan
params:
  name: 160-figureS1-ecology
subtitle: Redraw figS1_ecology from frozen records
title: 160-figureS1-ecology
---

This step reads the registered frozen inputs. It does not fit models,
run bootstrap resampling, call annotations or update the original
database. See the [output map](../metadata/result-code-map.tsv) and
[reproducibility scope](../docs/reproducibility.md).

## Tasks

``` r
pc057_draw_figure(step, "figS1_ecology")
```

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                   step status            output
    ## 1 160-figureS1-ecology   pass figS1_ecology.pdf
    ## 2 160-figureS1-ecology   pass figS1_ecology.png

<figure>
<img src="data/160-figureS1-ecology/figS1_ecology.png"
alt="plot of chunk figure" />
<figcaption aria-hidden="true">plot of chunk figure</figcaption>
</figure>
