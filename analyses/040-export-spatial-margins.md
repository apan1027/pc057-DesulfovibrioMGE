# 040-export-spatial-margins
Jinlong Ru, Cunli Pan

## Tasks

``` r
pc057_export_spatial_margins(step)
```

    ## [
    ##   {
    ##     "scope": "all_loci",
    ##     "distance_bp": 0,
    ##     "n_strata": 41,
    ##     "n_loci": 42,
    ##     "observed": 13,
    ##     "expected_from_margins": 10.514443509546936,
    ##     "expected_saved": 10.5144435095469,
    ##     "match": true
    ##   },
    ##   {
    ##     "scope": "all_loci",
    ##     "distance_bp": 5000,
    ##     "n_strata": 41,
    ##     "n_loci": 42,
    ##     "observed": 27,
    ##     "expected_from_margins": 13.02029626252406,
    ##     "expected_saved": 13.0202962625241,
    ##     "match": true
    ##   },
    ##   {
    ##     "scope": "all_loci",
    ##     "distance_bp": 10000,
    ##     "n_strata": 41,
    ##     "n_loci": 42,
    ##     "observed": 28,
    ##     "expected_from_margins": 15.756755726490864,
    ##     "expected_saved": 15.7567557264909,
    ##     "match": true
    ##   },
    ##   {
    ##     "scope": "clone_representatives",
    ##     "distance_bp": 0,
    ##     "n_strata": 18,
    ##     "n_loci": 19,
    ##     "observed": 6,
    ##     "expected_from_margins": 4.36331160436706,
    ##     "expected_saved": 4.36331160436706,
    ##     "match": true
    ##   },
    ##   {
    ##     "scope": "clone_representatives",
    ##     "distance_bp": 5000,
    ##     "n_strata": 18,
    ##     "n_loci": 19,
    ##     "observed": 11,
    ##     "expected_from_margins": 5.59193793423553,
    ##     "expected_saved": 5.59193793423553,
    ##     "match": true
    ##   },
    ##   {
    ##     "scope": "clone_representatives",
    ##     "distance_bp": 10000,
    ##     "n_strata": 18,
    ##     "n_loci": 19,
    ##     "observed": 12,
    ##     "expected_from_margins": 6.867726782438251,
    ##     "expected_saved": 6.86772678243825,
    ##     "match": true
    ##   }
    ## ]

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                         step status              output
    ## 1 040-export-spatial-margins   pass spatial-margins.csv
