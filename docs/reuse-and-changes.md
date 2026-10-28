# Code reuse and changes

Jinlong Ru authored the original PC057 R/Quarto project. The downloaded project metadata identifies `rujinlong/pc057-DesulfovibrioMGE` as its original GitHub source. Local snapshots preserve 29 numbered analysis QMD files, the original helper and accompanying scripts; this release does not infer the current remote repository state from that metadata.

The original MIT notice is preserved verbatim. The numbered notebook/path organisation follows the author's qproj framework. The Swedish-water project supplied an organisational reference; the recent RecWaterVirome project supplied examples of compact R/analyses/dev layout and English documentation. Their scientific data, conclusions, author lists and figure algorithms were not imported.

The current scientific revision separated Defense and Antidefense before graph merging, audited specimen metadata and used genome-conditioned and same-contig comparisons. Its source files were already saved before this release. This release parameterises their input/output locations and documents their relation to the manuscript. The exact source hashes and modifications are in `metadata/code-provenance.json`.

The old mixed-activity, host-only ecological and database-writing routes are not active entry points. `metadata/original-module-disposition.tsv` records which original modules were superseded, treated as annotation provenance or excluded as unrelated historical work. Exclusion from this compact release does not delete the source workspace or Git history of another repository.

Frozen plotting expressions are inherited from the validated saved-result redraw script. They were split into figure-specific functions and keep a fixed display seed. S2 is explicitly copied from its archived canonical image; its marker tree source is retained under `dev/upstream-source/`, with missing sequence/reference inputs documented.

ARG/VF integration rules are retained as standalone functions extracted from the original module. Only comments and function input/output plumbing were changed; database writes and unrelated historical statistical analyses were removed. They are not invoked by the default frozen workflow or claimed as newly validated annotation runs.

All new public-facing prose, comments and runtime messages are English. Original identifiers, specimen metadata and archived record values are preserved. Internal review conversations, private configuration and manuscript drafts are outside this release.
