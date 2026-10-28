# 020-export-tables
Jinlong Ru, Cunli Pan

## Tasks

``` r
pc057_export_tables(step)
```

    ## [
    ##   {
    ##     "table": "S1_Metadata",
    ##     "rows": 436,
    ##     "columns": 38,
    ##     "source": "metadata_corrected.csv"
    ##   },
    ##   {
    ##     "table": "S2_Load",
    ##     "rows": 436,
    ##     "columns": 20,
    ##     "source": "per_genome_corrected.csv, pivoted to one row per genome"
    ##   },
    ##   {
    ##     "table": "S3_DefenseCalls",
    ##     "rows": 6937,
    ##     "columns": 22,
    ##     "source": "calls_classified.csv, Defense + Defense_candidate"
    ##   },
    ##   {
    ##     "table": "S4_AntiDefense",
    ##     "rows": 493,
    ##     "columns": 22,
    ##     "source": "calls_classified.csv, Antidefense only"
    ##   },
    ##   {
    ##     "table": "S5_MergedDefense",
    ##     "rows": 5217,
    ##     "columns": 20,
    ##     "source": "call_tiers.csv, T3 defense"
    ##   },
    ##   {
    ##     "table": "S6_Enrichment",
    ##     "rows": 145,
    ##     "columns": 19,
    ##     "source": "enrichment_corrected.csv"
    ##   },
    ##   {
    ##     "table": "S7_ClusterIntervals",
    ##     "rows": 290,
    ##     "columns": 15,
    ##     "source": "cluster_intervals_corrected.csv"
    ##   },
    ##   {
    ##     "table": "S8_PDC_S31",
    ##     "rows": 42,
    ##     "columns": 24,
    ##     "source": "pdcs31_call_evidence.csv; raw PADLOC matches and specimen records"
    ##   },
    ##   {
    ##     "table": "S9_SpatialNull",
    ##     "rows": 6,
    ##     "columns": 8,
    ##     "source": "pdcs31_gene_position_null.csv"
    ##   },
    ##   {
    ##     "table": "S10_EcologyModels",
    ##     "rows": 33,
    ##     "columns": 15,
    ##     "source": "ecology_models_corrected.csv + model_diagnostics.csv; Wald intervals, nominal P; NB boundary flagged"
    ##   },
    ##   {
    ##     "table": "S11_Outgroup",
    ##     "rows": 1779,
    ##     "columns": 18,
    ##     "source": "outgroup_corrected.csv; 593 genomes x3 conventions; use qc_pass and ecology scope"
    ##   },
    ##   {
    ##     "table": "S12_OutgroupTests",
    ##     "rows": 48,
    ##     "columns": 8,
    ##     "source": "outgroup_tests_corrected.csv"
    ##   },
    ##   {
    ##     "table": "S13_OutgroupModels",
    ##     "rows": 18,
    ##     "columns": 12,
    ##     "source": "outgroup_models_corrected.csv; genome-level exploratory adjustment"
    ##   },
    ##   {
    ##     "table": "S14_ARG",
    ##     "rows": 1198,
    ##     "columns": 11,
    ##     "source": "archived arg_consensus, unchanged"
    ##   },
    ##   {
    ##     "table": "S15_VF_High",
    ##     "rows": 12782,
    ##     "columns": 33,
    ##     "source": "archived vf_consensus, HIGH tier, unchanged"
    ##   },
    ##   {
    ##     "table": "S16_ViralCatalog",
    ##     "rows": 3461,
    ##     "columns": 17,
    ##     "source": "archived viral_catalog + checkv_quality, unchanged"
    ##   },
    ##   {
    ##     "table": "S17_ViralLabels",
    ##     "rows": 3325,
    ##     "columns": 27,
    ##     "source": "archived vcontact3_assignments, unchanged"
    ##   },
    ##   {
    ##     "table": "S18_MetadataAudit",
    ##     "rows": 436,
    ##     "columns": 13,
    ##     "source": "metadata_audit.csv; original/revised labels and source evidence"
    ##   },
    ##   {
    ##     "table": "S19_SpatialStrata",
    ##     "rows": 177,
    ##     "columns": 7,
    ##     "source": "pdcs31_gene_position_strata.csv; full/clone sets at 0, 5 and 10 kb; 177 strata"
    ##   },
    ##   {
    ##     "table": "S20_SharedSpecies",
    ##     "rows": 53,
    ##     "columns": 5,
    ##     "source": "shared_species_corrected.csv; assigned species only"
    ##   },
    ##   {
    ##     "table": "S21_NeighborhoodGenes",
    ##     "rows": 100,
    ##     "columns": 17,
    ##     "source": "neighborhood_figure_genes.csv; Figure4 representative windows"
    ##   }
    ## ]

``` r
receipt <- pc057_step_finish(step)
print(receipt)
```

    ##                 step status                            output
    ## 1  020-export-tables   pass            tables/S1_Metadata.csv
    ## 2  020-export-tables   pass      tables/S10_EcologyModels.csv
    ## 3  020-export-tables   pass           tables/S11_Outgroup.csv
    ## 4  020-export-tables   pass      tables/S12_OutgroupTests.csv
    ## 5  020-export-tables   pass     tables/S13_OutgroupModels.csv
    ## 6  020-export-tables   pass                tables/S14_ARG.csv
    ## 7  020-export-tables   pass            tables/S15_VF_High.csv
    ## 8  020-export-tables   pass       tables/S16_ViralCatalog.csv
    ## 9  020-export-tables   pass        tables/S17_ViralLabels.csv
    ## 10 020-export-tables   pass      tables/S18_MetadataAudit.csv
    ## 11 020-export-tables   pass      tables/S19_SpatialStrata.csv
    ## 12 020-export-tables   pass                tables/S2_Load.csv
    ## 13 020-export-tables   pass      tables/S20_SharedSpecies.csv
    ## 14 020-export-tables   pass  tables/S21_NeighborhoodGenes.csv
    ## 15 020-export-tables   pass        tables/S3_DefenseCalls.csv
    ## 16 020-export-tables   pass         tables/S4_AntiDefense.csv
    ## 17 020-export-tables   pass       tables/S5_MergedDefense.csv
    ## 18 020-export-tables   pass          tables/S6_Enrichment.csv
    ## 19 020-export-tables   pass    tables/S7_ClusterIntervals.csv
    ## 20 020-export-tables   pass             tables/S8_PDC_S31.csv
    ## 21 020-export-tables   pass         tables/S9_SpatialNull.csv
    ## 22 020-export-tables   pass tables/table-export-manifest.json
