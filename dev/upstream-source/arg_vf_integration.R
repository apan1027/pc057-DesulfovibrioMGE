# ARG and VF integration expressions from Jinlong Ru's original module 06.
# Function interfaces replace global state; database writes are omitted.
# Annotation source records are external; these functions are not run by the default workflow.
# Requires dplyr, data.table, stringr and tibble.

pc057_integrate_arg <- function(df_abricate, df_deeparg) {
df_dpg_with_coords <- df_deeparg %>%
  dplyr::filter(!is.na(contig)) %>%
  dplyr::select(arg_id_dpg, ARG, contig, gene_start, gene_stop, locus_tag,
                ARO_name, Cut_Off, resistance_to_drug_classes_names)

dt_abr <- data.table::as.data.table(
  df_abricate[, c("arg_id_abr", "contig", "start", "stop", "gene_name", "resistance")]
)
data.table::setkey(dt_abr, contig, start, stop)

dt_dpg <- data.table::as.data.table(df_dpg_with_coords)
dt_dpg[, c("start", "stop") := .(gene_start, gene_stop)]

dt_overlap <- data.table::foverlaps(
  dt_dpg[, .(contig, start, stop, arg_id_dpg, ARG, ARO_name, Cut_Off,
             resistance_to_drug_classes_names)],
  dt_abr,
  by.x = c("contig", "start", "stop"),
  by.y = c("contig", "start", "stop"),
  type = "any", nomatch = NULL
)

# Build consensus table
df_consensus_both <- dt_overlap %>%
  as.data.frame() %>%
  dplyr::distinct(arg_id_abr, arg_id_dpg, .keep_all = TRUE) %>%
  dplyr::transmute(
    arg_id = paste0("CONS_", dplyr::row_number()),
    tools = "both",
    gene_name = paste0(gene_name, "/", ARG),
    drug_class = resistance_to_drug_classes_names,
    contig, start, stop = i.stop,
    arg_id_abr, arg_id_dpg,
    confidence_tier = "HIGH"
  )

# abricate-only
df_abr_only_ids <- setdiff(df_abricate$arg_id_abr, df_consensus_both$arg_id_abr)
df_abr_only <- df_abricate %>%
  dplyr::filter(arg_id_abr %in% df_abr_only_ids) %>%
  dplyr::transmute(
    arg_id = paste0("CONS_", nrow(df_consensus_both) + dplyr::row_number()),
    tools = "abricate",
    gene_name, drug_class = resistance,
    contig, start, stop,
    arg_id_abr, arg_id_dpg = NA_character_,
    confidence_tier = "MEDIUM"
  )

# deeparg-only (with coords)
df_dpg_only_ids <- setdiff(df_dpg_with_coords$arg_id_dpg, df_consensus_both$arg_id_dpg)
df_dpg_only <- df_dpg_with_coords %>%
  dplyr::filter(arg_id_dpg %in% df_dpg_only_ids) %>%
  dplyr::transmute(
    arg_id = paste0("CONS_", nrow(df_consensus_both) + nrow(df_abr_only) +
                    dplyr::row_number()),
    tools = "deeparg",
    gene_name = ARG, drug_class = resistance_to_drug_classes_names,
    contig, start = gene_start, stop = gene_stop,
    arg_id_abr = NA_character_, arg_id_dpg,
    confidence_tier = ifelse(Cut_Off == "Strict", "MEDIUM", "CANDIDATE")
  )

df_arg_consensus <- dplyr::bind_rows(df_consensus_both, df_abr_only, df_dpg_only) %>%
  dplyr::mutate(genome_id = sub("_.*", "", contig))

return(df_arg_consensus)
}

pc057_integrate_vf <- function(df_vfdbA_raw, df_vfdbB_raw, bakta_lookup_vf) {
get_best_hit <- function(df) {
  df %>%
    dplyr::group_by(query) %>%
    dplyr::slice_max(bits, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}
df_a_bh <- get_best_hit(df_vfdbA_raw)
df_b_bh <- get_best_hit(df_vfdbB_raw)

# Tier classification
classify_tier <- function(df, set_label) {
  df %>%
    dplyr::mutate(
      tier = dplyr::case_when(
        pident >= 50 & qcov >= 0.7 & tcov >= 0.7 & evalue <= 1e-30 ~ "HIGH",
        pident >= 30 & qcov >= 0.5 & tcov >= 0.5 & evalue <= 1e-10 ~ "MEDIUM",
        evalue <= 1e-5 ~ "CANDIDATE",
        TRUE ~ "UNRELIABLE"
      ),
      vfdb_set = set_label
    )
}
df_a_tier <- classify_tier(df_a_bh, "A")
df_b_tier <- classify_tier(df_b_bh, "B")

# Apply filtering rule: setA accepts HIGH+MED, setB HIGH only
df_vfdbA_filtered <- df_a_tier %>% dplyr::filter(tier %in% c("HIGH", "MEDIUM"))
df_vfdbB_filtered <- df_b_tier %>% dplyr::filter(tier == "HIGH")

# Parse theader for VF metadata via regex
parse_vf_header <- function(theader) {
  tibble::tibble(
    vfg_id = stringr::str_extract(theader, "^VFG\\d+"),
    gene_name = stringr::str_match(theader, "\\)\\s+\\(([^)]+)\\)\\s+")[, 2],
    vf_name = stringr::str_match(theader, "\\)\\s+\\([^)]+\\)\\s+([^\\[]+?)\\s*\\[")[, 2],
    vf_id = paste0("VF", stringr::str_match(theader, "\\(VF(\\d+)\\)")[, 2]),
    vfc_category = stringr::str_trim(
      stringr::str_match(theader, "\\)\\s*-\\s*([^()\\[]+?)\\s*\\(VFC")[, 2]),
    vfc_id = paste0("VFC", stringr::str_match(theader, "\\(VFC(\\d+)\\)")[, 2]),
    source_species = stringr::str_match(theader, "\\[([^\\[\\]]+)\\]\\s*$")[, 2]
  )
}

df_vfdbA_filtered <- df_vfdbA_filtered %>%
  dplyr::bind_cols(parse_vf_header(.$theader)) %>%
  dplyr::rename(locus_tag = query)
df_vfdbB_filtered <- df_vfdbB_filtered %>%
  dplyr::bind_cols(parse_vf_header(.$theader)) %>%
  dplyr::rename(locus_tag = query)

# Build consensus VF set: setA preferred; setB-only added
df_b_only <- df_vfdbB_filtered %>%
  dplyr::filter(!locus_tag %in% df_vfdbA_filtered$locus_tag)
df_vf_consensus <- dplyr::bind_rows(
  df_vfdbA_filtered %>% dplyr::mutate(source_set = "A"),
  df_b_only %>% dplyr::mutate(source_set = "B_only")
)


dplyr::left_join(df_vf_consensus,bakta_lookup_vf,by="locus_tag")
}
