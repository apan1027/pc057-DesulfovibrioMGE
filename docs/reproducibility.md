# Reproducibility levels

| Level | Included capability | Current validation |
| --- | --- | --- |
| Frozen result verification and exports | 72 input hashes; catalogue counts; ten RDS/CSV comparisons; all 21 table contents; six spatial-margin sets | Executed in this release |
| Figure display | Regenerate Figures 1–5 and S1 from saved data; copy canonical S2 | Executed; seven previews visually inspected |
| Reconstruction from annotation records | Parameterised read-only SQLite preparation, separated call-tier graph construction, and scientific-revision statistical source | Source parsed and transformations reviewed; not executed in this release |
| Original sequence-to-annotation pipeline | Original FASTA, caller environments, reference databases and external vpipe workflow | Not supplied or rerun; no end-to-end capability claimed |

Table-content validation used the manuscript's `Supplementary_Data_1_Tables.xlsx` as its reference, including headers. Its SHA-256 and the expected hashes of the 21 exported CSV files are recorded in [validation](../metadata/validation.json). The isolated-checkout comparison used the original executed checkout as its reference and the same existing R library. Later report-header and documentation changes did not rerun the full driver.

The optional `dev/scientific-source/recalculate.R` refuses to run without `--allow-refit` and a fresh output directory. It contains the original conditional tests, bootstrap and models, with only path/output plumbing and the already adopted S7/S19 reporting corrections changed. Invoking it starts inferential work and is separate from the tested default workflow. It is not executed by the report driver.

Optional examples, requiring the original annotation SQLite database:

```sh
python3 dev/scientific-source/prepare_data.py --database /path/to/pc057.sqlite --bundle analyses/data/00-raw/d00-resource/pc057 --output /path/to/new-prepared-records
python3 dev/scientific-source/build_call_tiers.py --input /path/to/new-prepared-records --output /path/to/new-tier-records
Rscript --vanilla dev/scientific-source/recalculate.R --input /path/to/new-tier-records --output /path/to/new-statistical-run --allow-refit
```

These are documented interfaces, not a claim that this sequence was validated in a clean environment. SQLite is opened read-only by preparation. The complete raw-sequence and annotation workflow is not reconstructed by these scripts.

S7's Defense carrier frames contain 436 genomes and 62 assigned species; Antidefense frames contain 271 and 46. Adding zero-call clusters would define another bootstrap design. S19 reports 177 strata across both full and clone-representative sets at 0, 5 and 10 kb. Statistical estimates are inherited, not recalculated during release preparation.

The lockfile fixes the observed R package graph and qproj commit. It does not lock Python, system libraries, external annotation databases or fonts. Python helpers require only Python's standard library. The test uses an existing R library; fresh-machine dependency restoration has not been verified. Native Quarto completion was blocked by temporary-directory cleanup restrictions; all eleven notebook R bodies and Markdown reports were executed with the supplied knitr/Pandoc driver.
