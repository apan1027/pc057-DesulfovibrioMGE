# 010-verify-inputs
Jinlong Ru, Cunli Pan

## Tasks

``` r
pc057_verify_inputs(step)
```

    ## {
    ##   "status": "pass",
    ##   "input_files": 72,
    ##   "core_genomes": 436,
    ##   "defense_calls": 6937,
    ##   "antidefense_calls": 493,
    ##   "T3_defense_candidates": 5217,
    ##   "spatial_strata": 177,
    ##   "spatial_margin_checks": [
    ##     {
    ##       "scope": "all_loci",
    ##       "distance_bp": 0,
    ##       "n": 42,
    ##       "observed": 13,
    ##       "expected": 10.514443509546936
    ##     },
    ##     {
    ##       "scope": "all_loci",
    ##       "distance_bp": 5000,
    ##       "n": 42,
    ##       "observed": 27,
    ##       "expected": 13.02029626252406
    ##     },
    ##     {
    ##       "scope": "all_loci",
    ##       "distance_bp": 10000,
    ##       "n": 42,
    ##       "observed": 28,
    ##       "expected": 15.756755726490864
    ##     },
    ##     {
    ##       "scope": "clone_representatives",
    ##       "distance_bp": 0,
    ##       "n": 19,
    ##       "observed": 6,
    ##       "expected": 4.36331160436706
    ##     },
    ##     {
    ##       "scope": "clone_representatives",
    ##       "distance_bp": 5000,
    ##       "n": 19,
    ##       "observed": 11,
    ##       "expected": 5.59193793423553
    ##     },
    ##     {
    ##       "scope": "clone_representatives",
    ##       "distance_bp": 10000,
    ##       "n": 19,
    ##       "observed": 12,
    ##       "expected": 6.867726782438251
    ##     }
    ##   ],
    ##   "no_inference": true
    ## }

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                step status                  output
    ## 1 010-verify-inputs   pass input-verification.json
