# Core data dictionary

| Field or concept | Meaning |
| --- | --- |
| `genome_id` | Internal core G-series or selected outgroup O-series record identifier; use the corresponding public assembly register |
| `accession` / `assembly_accession` | Public NCBI assembly accession, including its version |
| `gtdb_species` | Retained GTDB species assignment; blank values remain unassigned rather than a pooled species |
| `ecology_category` | Revised specimen-source category, supported by specimen fields and the metadata audit |
| `activity` | Original DefenseFinder Defense/Antidefense label or the explicitly putative PADLOC defense-candidate status |
| `category` | Separate Defense or Antidefense analysis category |
| `tier` in call tables | T2 deduplicated calls, T3 cross-tool overlap merged only for the same normalised name, or T4 any-overlap locus bound |
| `tier` in ARG/VF tables | Tool-specific annotation integration criteria; not the T2/T3/T4 convention or a shared calibrated phenotype-confidence scale |
| `near` / prophage context | Overlap with a predicted consensus prophage or separation of at most 10 kb |
| `context_class` | Fully contained, boundary overlap, external 0–5 kb, external 5–10 kb, or no context within 10 kb |
| `clone_cluster` | Retained provenance cluster used by the quality-first clone representative rule |
| `bioproject` | NCBI BioProject provenance; the ecological contrast is between studies in this collection |
| `N`, `K`, `n`, `x` in S19 | Eligible CDS positions, positions meeting the specified distance, observed candidate count and observed in-context candidate count for one contig |
| `q_conditional` | Previously calculated BH-adjusted genome-conditioned test result; no new test is performed by release verification |
| `activity_carrier_genomes` | Bootstrap frame of genomes carrying at least one call of the activity being analysed |
| `quality_eligible` | Saved outgroup completeness/contamination eligibility flag; 95 of 157 selected records qualify |

See `metadata/input-schema.tsv` for file schemas and row counts. Original source columns are retained for traceability; dictionary descriptions do not replace the manuscript's Methods or the original public database definitions.
