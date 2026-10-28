---
author:
- name: Jinlong Ru
- name: Cunli Pan
params:
  name: 030-verify-saved-results
subtitle: Compare ten saved RDS result components with CSV exports
title: 030-verify-saved-results
---

This step reads the registered frozen inputs. It does not fit models,
run bootstrap resampling, call annotations or update the original
database. See the [output map](../metadata/result-code-map.tsv) and
[reproducibility scope](../docs/reproducibility.md).

## Tasks

``` r
pc057_verify_saved_results(step)
```

    ## meta: 436 rows / 38 CSV columns agree
    ## per_genome: 2616 rows / 58 CSV columns agree
    ## enrichment: 145 rows / 19 CSV columns agree
    ## cluster_intervals: 290 rows / 15 CSV columns agree
    ## spatial: 6 rows / 8 CSV columns agree
    ## ecology_models: 33 rows / 12 CSV columns agree
    ## outgroup: 1779 rows / 18 CSV columns agree
    ## outgroup_tests: 48 rows / 8 CSV columns agree
    ## outgroup_models: 18 rows / 12 CSV columns agree
    ## shared_species: 53 rows / 5 CSV columns agree
    ## PASS: all ten saved-result components match their CSV exports; no inferential computation executed.

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                       step status                         output
    ## 1 030-verify-saved-results   pass saved-results-verification.txt
