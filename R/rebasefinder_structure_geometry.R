.dnmb_rebasefinder_structure_pair_rules <- function() {
  rules <- base::data.frame(
    pair_id = c(
      "amino_mtase_SAM-IV", "c5_mtase_PCQ-ENV", "mmei_mtase_I-IV",
      "hsdr_motor_WA-WB", "hsdr_motor_WB-MIII",
      "resiii_motor_WA-WB", "resiii_motor_WB-MIII", "pld_HKD-HKD"
    ),
    motif_a = c(
      "SAM", "N5C-PC", "MmeI-I", "P-loop", "HsdR-WB",
      "ResIII-WA", "ResIII-WB", "PLD-HKD"
    ),
    motif_b = c(
      "Amino-IV", "N5C-ENV", "Amino-IV", "HsdR-WB", "HsdR-MIII",
      "ResIII-WB", "ResIII-MIII", "PLD-HKD"
    ),
    family_pattern = c(
      "^Type (I|II|III)$|^Unclassified", "^Type II$|^Unclassified", "^Type II$",
      "^Type I$", "^Type I$", "^Type III$", "^Type III$", "^Type (II|IV)$"
    ),
    roles = c("M,RM", "M,RM", "M,RM", "R,RM", "R,RM", "R,RM", "R,RM", "R,RM"),
    max_min_ca = c(12, 12, 12, 20, 20, 20, 20, 20),
    max_centroid_ca = c(20, 20, 20, 30, 30, 30, 30, 30),
    stringsAsFactors = FALSE
  )
  base::rbind(
    rules,
    base::data.frame(
      pair_id = c("hsdr_nuclease_motor_PD-WA", "resiii_nuclease_motor_PD-WA"),
      motif_a = c("HsdR-PD", "ResIII-PD"),
      motif_b = c("P-loop", "ResIII-WA"),
      family_pattern = c("^Type I$", "^Type III$"),
      roles = c("R,RM", "R,RM"),
      max_min_ca = c(65, 65),
      max_centroid_ca = c(85, 85),
      stringsAsFactors = FALSE
    )
  )
}

.dnmb_rebasefinder_pdb_residue_code <- function(x) {
  codes <- c(
    ALA = "A", ARG = "R", ASN = "N", ASP = "D", CYS = "C",
    GLN = "Q", GLU = "E", GLY = "G", HIS = "H", ILE = "I",
    LEU = "L", LYS = "K", MET = "M", PHE = "F", PRO = "P",
    SER = "S", THR = "T", TRP = "W", TYR = "Y", VAL = "V",
    MSE = "M", SEC = "C", PYL = "K"
  )
  out <- base::unname(codes[base::toupper(base::trimws(base::as.character(x)))])
  out[base::is.na(out)] <- "X"
  out
}

.dnmb_rebasefinder_read_pdb_ca <- function(path) {
  empty <- base::data.frame(
    chain = character(), residue_key = character(), res_seq = integer(),
    insertion_code = character(), residue = character(), x = numeric(),
    y = numeric(), z = numeric(), confidence = numeric(),
    stringsAsFactors = FALSE
  )
  if (base::is.null(path) || base::is.na(path) || !base::file.exists(path) ||
      !base::grepl("[.]pdb$", path, ignore.case = TRUE)) return(empty)
  lines <- base::readLines(path, warn = FALSE)
  first_end <- base::which(base::grepl("^ENDMDL", lines))
  if (base::length(first_end)) lines <- lines[base::seq_len(first_end[[1]] - 1L)]
  atom <- lines[base::grepl("^ATOM  ", lines)]
  if (!base::length(atom)) return(empty)
  atom_name <- base::trimws(base::substr(atom, 13L, 16L))
  altloc <- base::substr(atom, 17L, 17L)
  atom <- atom[atom_name == "CA" & altloc %in% c(" ", "A", "")]
  if (!base::length(atom)) return(empty)
  num <- function(x) base::suppressWarnings(base::as.numeric(base::trimws(x)))
  chain <- base::trimws(base::substr(atom, 22L, 22L))
  chain[!base::nzchar(chain)] <- "_"
  res_seq <- base::suppressWarnings(base::as.integer(base::trimws(base::substr(atom, 23L, 26L))))
  insertion <- base::trimws(base::substr(atom, 27L, 27L))
  key <- base::paste(chain, res_seq, insertion, sep = ":")
  keep <- !base::duplicated(key) & !base::is.na(res_seq)
  confidence <- num(base::substr(atom[keep], 61L, 66L))
  finite_confidence <- confidence[base::is.finite(confidence)]
  if (base::length(finite_confidence) && base::max(finite_confidence) <= 1.5) {
    confidence <- confidence * 100
  }
  base::data.frame(
    chain = chain[keep], residue_key = key[keep], res_seq = res_seq[keep],
    insertion_code = insertion[keep],
    residue = .dnmb_rebasefinder_pdb_residue_code(base::substr(atom[keep], 18L, 20L)),
    x = num(base::substr(atom[keep], 31L, 38L)),
    y = num(base::substr(atom[keep], 39L, 46L)),
    z = num(base::substr(atom[keep], 47L, 54L)),
    confidence = confidence,
    stringsAsFactors = FALSE
  )
}

.dnmb_rebasefinder_structure_provenance <- function(path) {
  if (base::is.null(path) || base::is.na(path) || !base::file.exists(path)) return("unavailable")
  path_text <- base::tolower(base::as.character(path)[[1]])
  header <- base::readLines(path, n = 80L, warn = FALSE)
  if (base::grepl("promod3", path_text, fixed = TRUE) ||
      base::any(base::grepl("promod3|homology model", header, ignore.case = TRUE))) {
    return("promod3_homology")
  }
  if (base::grepl("alphafold|esmfold|query_structures|predicted", path_text, perl = TRUE) ||
      base::any(base::grepl("alphafold|esmfold|predicted model", header, ignore.case = TRUE))) {
    return("predicted")
  }
  "experimental_or_user"
}

.dnmb_rebasefinder_predicted_structure <- function(path) {
  .dnmb_rebasefinder_structure_provenance(path) %in% c("promod3_homology", "predicted")
}

.dnmb_rebasefinder_align_structure_chain <- function(query_sequence, ca) {
  query_sequence <- .dnmb_rebasefinder_normalize_protein(query_sequence)
  if (base::is.na(query_sequence) || !base::nrow(ca)) return(NULL)
  chains <- base::unique(ca$chain)
  candidates <- base::lapply(chains, function(chain) {
    chain_ca <- ca[ca$chain == chain, , drop = FALSE]
    structure_sequence <- base::paste0(chain_ca$residue, collapse = "")
    if (!base::nzchar(structure_sequence)) return(NULL)
    if (base::identical(query_sequence, structure_sequence)) {
      mapping <- base::data.frame(
        query_pos = base::seq_len(base::nchar(query_sequence)),
        ca_index = base::seq_len(base::nrow(chain_ca)),
        stringsAsFactors = FALSE
      )
      identity <- 1
    } else {
      alignment <- base::suppressWarnings(Biostrings::pairwiseAlignment(
        Biostrings::AAString(query_sequence),
        Biostrings::AAString(structure_sequence),
        type = "global"
      ))
      aligned_query <- base::strsplit(base::as.character(Biostrings::alignedPattern(alignment)), "", fixed = TRUE)[[1]]
      aligned_structure <- base::strsplit(base::as.character(Biostrings::alignedSubject(alignment)), "", fixed = TRUE)[[1]]
      query_pos <- 0L
      structure_pos <- 0L
      mapped_query <- integer()
      mapped_structure <- integer()
      identical_residue <- logical()
      for (i in base::seq_along(aligned_query)) {
        if (aligned_query[[i]] != "-") query_pos <- query_pos + 1L
        if (aligned_structure[[i]] != "-") structure_pos <- structure_pos + 1L
        if (aligned_query[[i]] != "-" && aligned_structure[[i]] != "-") {
          mapped_query <- c(mapped_query, query_pos)
          mapped_structure <- c(mapped_structure, structure_pos)
          identical_residue <- c(
            identical_residue,
            aligned_query[[i]] == aligned_structure[[i]] ||
              aligned_query[[i]] == "X" || aligned_structure[[i]] == "X"
          )
        }
      }
      mapping <- base::data.frame(
        query_pos = mapped_query,
        ca_index = mapped_structure,
        stringsAsFactors = FALSE
      )
      identity <- if (base::length(identical_residue)) base::mean(identical_residue) else 0
    }
    coverage <- base::nrow(mapping) / base::nchar(query_sequence)
    list(
      chain = chain, ca = chain_ca, mapping = mapping,
      sequence_identity = identity, sequence_coverage = coverage,
      rank = coverage * identity
    )
  })
  candidates <- base::Filter(base::Negate(base::is.null), candidates)
  if (!base::length(candidates)) return(NULL)
  ranks <- base::vapply(candidates, `[[`, numeric(1), "rank")
  candidates[[base::which.max(ranks)]]
}

.dnmb_rebasefinder_motif_coordinates <- function(hit, alignment) {
  start <- base::suppressWarnings(base::as.integer(hit$start_aa[[1]]))
  end <- base::suppressWarnings(base::as.integer(hit$end_aa[[1]]))
  if (base::is.na(start) || base::is.na(end) || end < start) return(NULL)
  positions <- base::seq.int(start, end)
  map_idx <- base::match(positions, alignment$mapping$query_pos)
  ca_idx <- alignment$mapping$ca_index[map_idx[!base::is.na(map_idx)]]
  ca_idx <- ca_idx[!base::is.na(ca_idx) & ca_idx >= 1L & ca_idx <= base::nrow(alignment$ca)]
  coords <- alignment$ca[base::unique(ca_idx), , drop = FALSE]
  list(
    coords = coords,
    coverage = base::nrow(coords) / base::length(positions),
    mean_confidence = if (base::nrow(coords) && base::any(base::is.finite(coords$confidence))) {
      base::mean(coords$confidence[base::is.finite(coords$confidence)])
    } else {
      NA_real_
    }
  )
}

.dnmb_rebasefinder_ca_pair_geometry <- function(a, b) {
  if (base::is.null(a) || base::is.null(b) || !base::nrow(a$coords) || !base::nrow(b$coords)) {
    return(list(min_ca = NA_real_, centroid = NA_real_, contacts_8a = NA_integer_))
  }
  dx <- base::outer(a$coords$x, b$coords$x, "-")
  dy <- base::outer(a$coords$y, b$coords$y, "-")
  dz <- base::outer(a$coords$z, b$coords$z, "-")
  distances <- base::sqrt(dx^2 + dy^2 + dz^2)
  centroid_a <- base::colMeans(a$coords[, c("x", "y", "z"), drop = FALSE])
  centroid_b <- base::colMeans(b$coords[, c("x", "y", "z"), drop = FALSE])
  list(
    min_ca = base::min(distances, na.rm = TRUE),
    centroid = base::sqrt(base::sum((centroid_a - centroid_b)^2)),
    contacts_8a = base::sum(distances <= 8, na.rm = TRUE)
  )
}

.dnmb_rebasefinder_geometry_foldseek_state <- function(row) {
  supported <- if ("REBASEfinder_structure_supported" %in% base::names(row)) {
    row$REBASEfinder_structure_supported[[1]] %in% TRUE
  } else if ("REBASEfinder_structure_pass" %in% base::names(row)) {
    row$REBASEfinder_structure_pass[[1]] %in% TRUE
  } else {
    FALSE
  }
  consistent <- if ("REBASEfinder_structure_candidate_consistent" %in% base::names(row)) {
    row$REBASEfinder_structure_candidate_consistent[[1]]
  } else {
    NA
  }
  status <- if ("REBASEfinder_structure_status" %in% base::names(row)) {
    base::as.character(row$REBASEfinder_structure_status[[1]])
  } else {
    NA_character_
  }
  if (supported && !base::identical(consistent, FALSE)) return("supported")
  if (base::identical(consistent, FALSE) ||
      (!base::is.na(status) && base::grepl("mismatch|conflict", status, ignore.case = TRUE))) {
    return("conflict")
  }
  if (!base::is.na(status) && base::nzchar(status) && status != "structure_missing") return("weak")
  "unavailable"
}

.dnmb_rebasefinder_empty_geometry_table <- function() {
  base::data.frame(
    locus_tag = character(), family_id = character(), enzyme_role = character(),
    pair_id = character(), motif_a = character(), motif_b = character(),
    motif_a_start = integer(), motif_a_end = integer(), motif_b_start = integer(),
    motif_b_end = integer(), structure_path = character(), structure_chain = character(),
    structure_method = character(), predicted_structure = logical(), mapping_identity = numeric(), mapping_coverage = numeric(),
    motif_a_coverage = numeric(), motif_b_coverage = numeric(), motif_a_mean_plddt = numeric(),
    motif_b_mean_plddt = numeric(), min_ca_distance = numeric(), centroid_distance = numeric(),
    contact_count_8a = integer(), sequence_separation = integer(), max_min_ca = numeric(),
    max_centroid_ca = numeric(), geometry_status = character(), foldseek_state = character(),
    combined_status = character(), stringsAsFactors = FALSE
  )
}

.dnmb_rebasefinder_verify_motif_geometry <- function(motif_hits,
                                                     protein_table,
                                                     structure_dirs,
                                                     structure_paths = NULL,
                                                     min_mapping_coverage = 0.3,
                                                     min_mapping_identity = 0.7,
                                                     min_motif_coverage = 0.8,
                                                     min_plddt = 70) {
  motif_hits <- base::as.data.frame(motif_hits, stringsAsFactors = FALSE)
  protein_table <- base::as.data.frame(protein_table, stringsAsFactors = FALSE)
  if (!base::nrow(motif_hits) || !base::nrow(protein_table) ||
      !base::all(c("locus_tag", "motif", "start_aa", "end_aa") %in% base::names(motif_hits)) ||
      !base::all(c("locus_tag", "translation") %in% base::names(protein_table))) {
    return(list(pairs = .dnmb_rebasefinder_empty_geometry_table(), summary = base::data.frame()))
  }
  rules <- .dnmb_rebasefinder_structure_pair_rules()
  motif_hits <- motif_hits[!base::is.na(motif_hits$start_aa) & !base::is.na(motif_hits$end_aa), , drop = FALSE]
  if ("evidence_level" %in% base::names(motif_hits)) {
    motif_hits <- motif_hits[motif_hits$evidence_level %in% c("supported", "sequence_hint"), , drop = FALSE]
  }
  if (!base::nrow(motif_hits)) {
    return(list(pairs = .dnmb_rebasefinder_empty_geometry_table(), summary = base::data.frame()))
  }
  protein_table$locus_tag <- .dnmb_module_clean_annotation_key(protein_table$locus_tag)
  loci <- base::unique(.dnmb_module_clean_annotation_key(motif_hits$locus_tag))
  pair_rows <- list()
  n_pair <- 0L
  for (locus in loci) {
    hits <- motif_hits[.dnmb_module_clean_annotation_key(motif_hits$locus_tag) == locus, , drop = FALSE]
    protein_idx <- base::match(locus, protein_table$locus_tag)
    if (base::is.na(protein_idx)) next
    protein <- protein_table[protein_idx, , drop = FALSE]
    family <- if ("family_id" %in% base::names(hits)) base::as.character(hits$family_id[[1]]) else NA_character_
    rule_family <- .dnmb_rebasefinder_canonical_structure_family(family)
    role <- if ("enzyme_role" %in% base::names(hits)) base::as.character(hits$enzyme_role[[1]]) else NA_character_
    applicable <- rules[
      base::vapply(base::seq_len(base::nrow(rules)), function(i) {
        family_ok <- base::is.na(rule_family) || !base::nzchar(rule_family) ||
          base::grepl(rules$family_pattern[[i]], rule_family, ignore.case = TRUE)
        role_ok <- base::is.na(role) || !base::nzchar(role) ||
          role %in% base::trimws(base::strsplit(rules$roles[[i]], ",", fixed = TRUE)[[1]])
        family_ok && role_ok
      }, logical(1)), , drop = FALSE
    ]
    if (!base::nrow(applicable)) next
    structure_path <- .dnmb_rebasefinder_structure_file_for_query(
      locus,
      structure_dirs,
      explicit_paths = structure_paths
    )
    ca <- .dnmb_rebasefinder_read_pdb_ca(structure_path)
    alignment <- if (base::nrow(ca)) {
      .dnmb_rebasefinder_align_structure_chain(protein$translation[[1]], ca)
    } else {
      NULL
    }
    structure_method <- .dnmb_rebasefinder_structure_provenance(structure_path)
    predicted <- structure_method %in% c("promod3_homology", "predicted")
    foldseek_state <- .dnmb_rebasefinder_geometry_foldseek_state(protein)
    for (j in base::seq_len(base::nrow(applicable))) {
      rule <- applicable[j, , drop = FALSE]
      a_idx <- base::which(hits$motif == rule$motif_a[[1]])
      b_idx <- base::which(hits$motif == rule$motif_b[[1]])
      if (!base::length(a_idx) || !base::length(b_idx)) next
      combinations <- base::expand.grid(a = a_idx, b = b_idx, KEEP.OUT.ATTRS = FALSE)
      if (rule$motif_a[[1]] == rule$motif_b[[1]]) combinations <- combinations[combinations$a < combinations$b, , drop = FALSE]
      if (!base::nrow(combinations)) next
      candidate_rows <- list()
      for (k in base::seq_len(base::nrow(combinations))) {
        hit_a <- hits[combinations$a[[k]], , drop = FALSE]
        hit_b <- hits[combinations$b[[k]], , drop = FALSE]
        coord_a <- if (!base::is.null(alignment)) .dnmb_rebasefinder_motif_coordinates(hit_a, alignment) else NULL
        coord_b <- if (!base::is.null(alignment)) .dnmb_rebasefinder_motif_coordinates(hit_b, alignment) else NULL
        geometry <- .dnmb_rebasefinder_ca_pair_geometry(coord_a, coord_b)
        mapping_ok <- !base::is.null(alignment) &&
          alignment$sequence_coverage >= min_mapping_coverage &&
          alignment$sequence_identity >= min_mapping_identity
        coverage_ok <- !base::is.null(coord_a) && !base::is.null(coord_b) &&
          coord_a$coverage >= min_motif_coverage && coord_b$coverage >= min_motif_coverage
        confidence_values <- if (!base::is.null(coord_a) && !base::is.null(coord_b) &&
                                 structure_method != "promod3_homology") {
          c(coord_a$mean_confidence, coord_b$mean_confidence)
        } else {
          c(NA_real_, NA_real_)
        }
        confidence_ok <- structure_method == "promod3_homology" || !predicted ||
          (base::all(base::is.finite(confidence_values)) && base::min(confidence_values) >= min_plddt)
        geometry_status <- if (base::is.na(structure_path) || !base::file.exists(structure_path)) {
          "no_structure"
        } else if (!base::grepl("[.]pdb$", structure_path, ignore.case = TRUE)) {
          "unsupported_structure_format"
        } else if (!mapping_ok) {
          "mapping_failed"
        } else if (!coverage_ok) {
          "partial_motif_coverage"
        } else if (!confidence_ok) {
          "low_local_confidence"
        } else if (!base::is.na(geometry$min_ca) && !base::is.na(geometry$centroid) &&
                   geometry$min_ca <= rule$max_min_ca[[1]] &&
                   geometry$centroid <= rule$max_centroid_ca[[1]]) {
          "geometry_compatible"
        } else {
          "geometry_far"
        }
        combined_status <- if (geometry_status == "geometry_compatible" && structure_method == "promod3_homology") {
          "homology_model_supported"
        } else if (geometry_status == "geometry_compatible" && foldseek_state == "supported") {
          "3d_fold_supported"
        } else if (geometry_status == "geometry_compatible" && foldseek_state == "conflict") {
          "3d_fold_conflict"
        } else if (geometry_status == "geometry_compatible") {
          "3d_geometry_supported"
        } else {
          geometry_status
        }
        candidate_rows[[k]] <- base::data.frame(
          locus_tag = locus, family_id = family, enzyme_role = role,
          pair_id = rule$pair_id[[1]], motif_a = rule$motif_a[[1]], motif_b = rule$motif_b[[1]],
          motif_a_start = base::as.integer(hit_a$start_aa[[1]]), motif_a_end = base::as.integer(hit_a$end_aa[[1]]),
          motif_b_start = base::as.integer(hit_b$start_aa[[1]]), motif_b_end = base::as.integer(hit_b$end_aa[[1]]),
          structure_path = structure_path,
          structure_chain = if (!base::is.null(alignment)) alignment$chain else NA_character_,
          structure_method = structure_method,
          predicted_structure = predicted,
          mapping_identity = if (!base::is.null(alignment)) alignment$sequence_identity else NA_real_,
          mapping_coverage = if (!base::is.null(alignment)) alignment$sequence_coverage else NA_real_,
          motif_a_coverage = if (!base::is.null(coord_a)) coord_a$coverage else NA_real_,
          motif_b_coverage = if (!base::is.null(coord_b)) coord_b$coverage else NA_real_,
          motif_a_mean_plddt = if (structure_method != "promod3_homology" && !base::is.null(coord_a)) coord_a$mean_confidence else NA_real_,
          motif_b_mean_plddt = if (structure_method != "promod3_homology" && !base::is.null(coord_b)) coord_b$mean_confidence else NA_real_,
          min_ca_distance = geometry$min_ca, centroid_distance = geometry$centroid,
          contact_count_8a = geometry$contacts_8a,
          sequence_separation = base::abs(base::as.integer(hit_a$start_aa[[1]]) - base::as.integer(hit_b$start_aa[[1]])),
          max_min_ca = rule$max_min_ca[[1]], max_centroid_ca = rule$max_centroid_ca[[1]],
          geometry_status = geometry_status, foldseek_state = foldseek_state,
          combined_status = combined_status,
          stringsAsFactors = FALSE
        )
      }
      candidate_tbl <- base::do.call(base::rbind, candidate_rows)
      rank <- base::match(candidate_tbl$combined_status, c(
        "homology_model_supported", "3d_fold_supported", "3d_geometry_supported", "3d_fold_conflict",
        "low_local_confidence", "partial_motif_coverage", "geometry_far",
        "mapping_failed", "unsupported_structure_format", "no_structure"
      ))
      rank[base::is.na(rank)] <- 99L
      candidate_tbl <- candidate_tbl[base::order(
        rank,
        base::ifelse(base::is.na(candidate_tbl$min_ca_distance), Inf, candidate_tbl$min_ca_distance)
      ), , drop = FALSE]
      n_pair <- n_pair + 1L
      pair_rows[[n_pair]] <- candidate_tbl[1, , drop = FALSE]
    }
  }
  pairs <- if (base::length(pair_rows)) base::do.call(base::rbind, pair_rows) else .dnmb_rebasefinder_empty_geometry_table()
  if (!base::nrow(pairs)) return(list(pairs = pairs, summary = base::data.frame()))
  summary_rows <- base::lapply(base::split(pairs, pairs$locus_tag), function(x) {
    compatible <- x$combined_status %in% c(
      "homology_model_supported", "3d_fold_supported", "3d_geometry_supported"
    )
    n_compatible <- base::sum(compatible)
    homology_model <- base::any(x$structure_method == "promod3_homology")
    status <- if (homology_model && base::all(x$combined_status == "homology_model_supported")) {
      "homology_model_supported"
    } else if (homology_model && base::any(compatible)) {
      "homology_model_partial"
    } else if (base::all(x$combined_status == "3d_fold_supported")) {
      "3d_fold_supported"
    } else if (base::all(compatible)) {
      "3d_geometry_supported"
    } else if (base::any(compatible)) {
      "3d_partial"
    } else if (base::any(x$combined_status == "3d_fold_conflict")) {
      "3d_fold_conflict"
    } else if (base::any(x$combined_status == "geometry_far")) {
      if (homology_model) "homology_model_geometry_far" else "3d_far"
    } else if (base::any(x$combined_status == "low_local_confidence")) {
      "low_local_confidence"
    } else if (base::any(x$combined_status == "mapping_failed")) {
      "mapping_failed"
    } else if (base::all(x$combined_status == "no_structure")) {
      "no_structure"
    } else {
      x$combined_status[[1]]
    }
    canonical_family <- .dnmb_rebasefinder_canonical_structure_family(x$family_id[[1]])
    restriction_motor <- x$enzyme_role[[1]] %in% c("R", "RM") &&
      canonical_family %in% c("Type I", "Type III")
    nuclease_motor_pair <- base::any(x$pair_id %in%
      c("hsdr_nuclease_motor_PD-WA", "resiii_nuclease_motor_PD-WA"))
    if (restriction_motor && !nuclease_motor_pair && status %in%
        c("homology_model_supported", "homology_model_partial", "3d_fold_supported", "3d_geometry_supported", "3d_partial")) {
      status <- if (homology_model) "homology_model_motor_only" else "3d_motor_only"
    }
    mapping_coverage <- if (base::any(base::is.finite(x$mapping_coverage))) {
      base::max(x$mapping_coverage, na.rm = TRUE)
    } else {
      NA_real_
    }
    model_scope <- if (!base::is.finite(mapping_coverage)) {
      "coverage_unknown"
    } else if (mapping_coverage >= 0.80) {
      "near_full_model"
    } else {
      "partial_model"
    }
    label <- base::switch(
      status,
      `homology_model_supported` = if (model_scope == "partial_model") {
        base::sprintf("HOM\nlocal\n%d/%d", n_compatible, base::nrow(x))
      } else {
        base::sprintf("HOM\n%d/%d", n_compatible, base::nrow(x))
      },
      `homology_model_partial` = base::sprintf("HOM partial\n%d/%d", n_compatible, base::nrow(x)),
      `homology_model_motor_only` = "HOM motor\nPD?",
      `homology_model_geometry_far` = "HOM far",
      `3d_fold_supported` = if (model_scope == "partial_model") {
        base::sprintf("3D+FS\nlocal\n%d/%d", n_compatible, base::nrow(x))
      } else {
        base::sprintf("3D+FS\n%d/%d", n_compatible, base::nrow(x))
      },
      `3d_geometry_supported` = if (model_scope == "partial_model") {
        base::sprintf("3D\nlocal\n%d/%d", n_compatible, base::nrow(x))
      } else {
        base::sprintf("3D\n%d/%d", n_compatible, base::nrow(x))
      },
      `3d_partial` = base::sprintf("3D partial\n%d/%d", n_compatible, base::nrow(x)),
      `3d_motor_only` = "3D motor\nPD?",
      `3d_fold_conflict` = "3D/FS\nconflict",
      `3d_far` = "3D far",
      `low_local_confidence` = "3D lowQ",
      `mapping_failed` = "3D map?",
      ""
    )
    base::data.frame(
      locus_tag = x$locus_tag[[1]], structural_adjacency_status = status,
      structural_adjacency_label = label, n_pairs = base::nrow(x),
      n_compatible = n_compatible,
      min_ca_distance = if (base::any(base::is.finite(x$min_ca_distance))) base::min(x$min_ca_distance, na.rm = TRUE) else NA_real_,
      structure_mapping_coverage = mapping_coverage,
      structure_model_scope = model_scope,
      structure_method = base::paste(base::unique(x$structure_method), collapse = ","),
      foldseek_state = base::paste(base::unique(x$foldseek_state), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
  summary <- base::do.call(base::rbind, summary_rows)
  base::rownames(pairs) <- NULL
  base::rownames(summary) <- NULL
  list(pairs = pairs, summary = summary)
}

.dnmb_rebasefinder_write_motif_geometry <- function(output_dir, tbl, motif_hits) {
  structure_path <- .dnmb_rebasefinder_structure_path(output_dir)
  structure_dirs <- .dnmb_rebasefinder_structure_dirs(output_dir, structure_path)
  explicit_paths <- .dnmb_rebasefinder_homology_structure_path_map(tbl)
  result <- .dnmb_rebasefinder_verify_motif_geometry(
    motif_hits = motif_hits,
    protein_table = tbl,
    structure_dirs = structure_dirs,
    structure_paths = explicit_paths
  )
  module_dir <- if (base::basename(base::normalizePath(output_dir, winslash = "/", mustWork = FALSE)) == "dnmb_module_rebasefinder") {
    output_dir
  } else {
    base::file.path(output_dir, "dnmb_module_rebasefinder")
  }
  base::dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
  tsv <- base::file.path(module_dir, "DNMB_REBASEfinder_motif_geometry.tsv")
  utils::write.table(result$pairs, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  xlsx <- base::file.path(module_dir, "DNMB_REBASEfinder_motif_geometry.xlsx")
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch(
      openxlsx::write.xlsx(
        list(Pair_geometry = result$pairs, Gene_summary = result$summary),
        xlsx,
        overwrite = TRUE
      ),
      error = function(e) NULL
    )
  }
  result$tsv <- base::normalizePath(tsv, winslash = "/", mustWork = FALSE)
  result$xlsx <- if (base::file.exists(xlsx)) base::normalizePath(xlsx, winslash = "/", mustWork = FALSE) else NA_character_
  result
}

.dnmb_rebasefinder_motif_contact_path <- function(output_dir) {
  output_dir <- base::normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  module_dir <- if (base::basename(output_dir) == "dnmb_module_rebasefinder") {
    output_dir
  } else {
    base::file.path(output_dir, "dnmb_module_rebasefinder")
  }
  candidates <- c(
    base::file.path(module_dir, "DNMB_REBASEfinder_motif_contacts.tsv"),
    base::file.path(output_dir, "DNMB_REBASEfinder_motif_contacts.tsv")
  )
  existing <- candidates[base::file.exists(candidates)]
  if (base::length(existing)) existing[[1]] else candidates[[1]]
}

.dnmb_rebasefinder_read_motif_contacts <- function(output_dir, write_xlsx = TRUE) {
  path <- .dnmb_rebasefinder_motif_contact_path(output_dir)
  empty <- list(pairs = base::data.frame(), summary = base::data.frame(), tsv = path, xlsx = NA_character_)
  if (!base::file.exists(path) || base::file.info(path)$size <= 0) return(empty)
  pairs <- tryCatch(
    utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = ""),
    error = function(e) NULL
  )
  required <- c("locus_tag", "pair_id", "contact_status")
  if (base::is.null(pairs) || !base::nrow(pairs) || !base::all(required %in% base::names(pairs))) {
    return(empty)
  }
  numeric_cols <- base::intersect(c(
    "max_contact_probability", "mean_contact_probability",
    "top3_mean_contact_probability", "separation_matched_percentile",
    "sequence_separation", "inference_seconds"
  ), base::names(pairs))
  for (col in numeric_cols) pairs[[col]] <- base::suppressWarnings(base::as.numeric(pairs[[col]]))
  pairs$locus_tag <- .dnmb_module_clean_annotation_key(pairs$locus_tag)
  pairs <- pairs[!base::is.na(pairs$locus_tag) & base::nzchar(pairs$locus_tag), , drop = FALSE]
  if (!base::nrow(pairs)) return(empty)

  summary_rows <- base::lapply(base::split(pairs, pairs$locus_tag), function(x) {
    supported <- x$contact_status %in% c("contact_strong", "contact_supportive")
    n_supported <- base::sum(supported)
    status <- if (base::any(x$contact_status == "contact_strong")) {
      "contact_strong"
    } else if (base::any(x$contact_status == "contact_supportive")) {
      "contact_supportive"
    } else if (base::all(x$contact_status == "contact_weak")) {
      "contact_weak"
    } else {
      "contact_unavailable"
    }
    label <- if (n_supported > 0L) {
      base::sprintf("CMAP+\n%d/%d", n_supported, base::nrow(x))
    } else {
      ""
    }
    max_probability <- if ("max_contact_probability" %in% base::names(x) &&
                           base::any(base::is.finite(x$max_contact_probability))) {
      base::max(x$max_contact_probability, na.rm = TRUE)
    } else {
      NA_real_
    }
    max_percentile <- if ("separation_matched_percentile" %in% base::names(x) &&
                          base::any(base::is.finite(x$separation_matched_percentile))) {
      base::max(x$separation_matched_percentile, na.rm = TRUE)
    } else {
      NA_real_
    }
    base::data.frame(
      locus_tag = x$locus_tag[[1]],
      motif_contact_status = status,
      motif_contact_label = label,
      n_pairs = base::nrow(x),
      n_supported = n_supported,
      max_contact_probability = max_probability,
      max_separation_matched_percentile = max_percentile,
      model = if ("model" %in% base::names(x)) base::paste(base::unique(x$model), collapse = ",") else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  summary <- base::do.call(base::rbind, summary_rows)
  base::rownames(summary) <- NULL
  xlsx <- base::sub("[.]tsv$", ".xlsx", path, ignore.case = TRUE)
  if (base::isTRUE(write_xlsx) && requireNamespace("openxlsx", quietly = TRUE)) {
    tryCatch(
      openxlsx::write.xlsx(list(Pair_contacts = pairs, Gene_summary = summary), xlsx, overwrite = TRUE),
      error = function(e) NULL
    )
  }
  list(
    pairs = pairs,
    summary = summary,
    tsv = base::normalizePath(path, winslash = "/", mustWork = FALSE),
    xlsx = if (base::file.exists(xlsx)) base::normalizePath(xlsx, winslash = "/", mustWork = FALSE) else NA_character_
  )
}
