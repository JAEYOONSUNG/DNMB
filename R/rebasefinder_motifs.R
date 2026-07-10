.dnmb_rebasefinder_normalize_protein <- function(sequence) {
  sequence <- base::as.character(sequence)[1]
  if (base::is.na(sequence) || !base::nzchar(sequence)) return(NA_character_)
  sequence <- base::toupper(base::gsub("[^A-Za-z]", "", sequence))
  if (!base::nzchar(sequence)) NA_character_ else sequence
}

.dnmb_rebasefinder_motif_catalog <- function() {
  list(
    # DNA methyltransferase motifs. Short motifs become supporting evidence only
    # after a coherent motif set and a DNA-MTase annotation/domain are present.
    "SAM" = list(
      pattern = "[FYW].G.[GA]", expected_role = "M", axis = "M",
      full = "AdoMet-binding motif I FxGxG/FxGxA",
      pos_range = NULL, weight = 3L, evidence_kind = "sequence"
    ),
    "Amino-IV" = list(
      pattern = "[DNS]PP[YFW]", expected_role = "M", axis = "M",
      full = "Amino-DNA-MTase catalytic motif IV [DNS]PP[YFW]",
      pos_range = NULL, weight = 3L, evidence_kind = "sequence"
    ),
    "N5C-PC" = list(
      pattern = "PCQ", expected_role = "M", axis = "M",
      full = "C5-DNA-MTase catalytic motif IV PCQ",
      pos_range = NULL, weight = 3L, evidence_kind = "sequence"
    ),
    "N5C-ENV" = list(
      pattern = "ENV", expected_role = "M", axis = "M",
      full = "C5-DNA-MTase proton-donor motif VI ENV",
      pos_range = NULL, weight = 2L, evidence_kind = "sequence"
    ),
    "N5C-QxRxR" = list(
      pattern = "Q.R.R", expected_role = "M", axis = "M",
      full = "C5-DNA-MTase motif VIII QxRxR",
      pos_range = NULL, weight = 1L, evidence_kind = "sequence"
    ),

    # MmeI/Type IIL enzymes use a gamma-class MTase motif order and an
    # N-terminal D-X(8-14)-E-X(0-3)-K nuclease center.
    "MmeI-X" = list(
      pattern = "GAHY[TS]", expected_role = "RM", expected_family = "Type II",
      axis = "M", full = "MmeI-family gamma-class MTase motif X GAHY[TS]",
      pos_range = c(0.03, 0.55), weight = 2L, evidence_kind = "sequence"
    ),
    "MmeI-I" = list(
      pattern = "F[FYLIV]DPACG[CS]GNF", expected_role = "RM", expected_family = "Type II",
      axis = "M", full = "MmeI-family AdoMet-binding motif I F[FLIVY]DPACG[CS]GNF",
      pos_range = c(0.05, 0.65), weight = 3L, evidence_kind = "sequence"
    ),
    "MmeI-PDExK" = list(
      pattern = "D.{8,14}E.{0,3}K", expected_role = "RM", expected_family = "Type II",
      axis = "R", full = "MmeI-family N-terminal nuclease D-X(8-14)-E-X(0-3)-K",
      pos_range = c(0.00, 0.30), weight = 3L, evidence_kind = "sequence"
    ),
    "MmeI-architecture" = list(
      pattern = NA_character_, expected_role = "RM", expected_family = "Type II",
      axis = "RM", full = "MmeI-family fused REase-MTase-TRD domain architecture",
      pos_range = NULL, weight = 3L, evidence_kind = "domain_architecture"
    ),

    # Type I HsdR: N-terminal nuclease followed by an SF2 motor. The short
    # motor patterns never verify a protein independently.
    "HsdR-PD" = list(
      pattern = "PD.{5,30}[DE].K", expected_role = "R", expected_family = "Type I",
      axis = "R", full = "Type I HsdR N-terminal PD-(D/E)xK-like nuclease motif",
      pos_range = c(0.00, 0.40), weight = 3L, evidence_kind = "sequence"
    ),
    "P-loop" = list(
      pattern = "(?:[AG].{4}GK[ST]|[GAK]SGK[ST]|AAGK[ST])", expected_role = "R", expected_family = "Type I",
      axis = "R", full = "Type I HsdR Walker-A/P-loop motif",
      pos_range = c(0.10, 0.80), weight = 2L, evidence_kind = "sequence"
    ),
    "HsdR-WB" = list(
      pattern = "DE[AHCFQ][DHQ]", expected_role = "R", expected_family = "Type I",
      axis = "R", full = "Type I HsdR Walker-B/motif-II DEAD- or DExH-like sequence",
      pos_range = c(0.10, 0.85), weight = 2L, evidence_kind = "sequence"
    ),
    "HsdR-MIII" = list(
      pattern = "[ST]AT", expected_role = "R", expected_family = "Type I",
      axis = "R", full = "Type I HsdR SF2 motor motif III [ST]AT",
      pos_range = c(0.10, 0.85), weight = 1L, evidence_kind = "sequence"
    ),
    "HsdS-2TRD" = list(
      pattern = NA_character_, expected_role = "S", expected_family = "Type I",
      axis = "S", full = "Type I HsdS dual target-recognition-domain architecture",
      pos_range = NULL, weight = 3L, evidence_kind = "domain_architecture"
    ),

    # Alternative Type II nuclease folds. HNH and PLD/HKD are accepted only
    # with their matching domain annotations; one HKD motif is insufficient.
    "PD-ExK" = list(
      pattern = "PD.{5,30}[DE].K", expected_role = "R", expected_family = "Type II",
      axis = "R", full = "Type II REase PD-(D/E)xK-like catalytic motif",
      pos_range = NULL, weight = 3L, evidence_kind = "sequence"
    ),
    "HNH" = list(
      pattern = "H.{1,3}N.{5,25}H", expected_role = "R", expected_family = "Type II",
      axis = "R", full = "HNH-like nuclease sequence candidate (domain-gated)",
      pos_range = NULL, weight = 2L, evidence_kind = "sequence"
    ),
    "GIY-YIG" = list(
      pattern = "GIY.{5,20}YIG", expected_role = "R", expected_family = "Type II",
      axis = "R", full = "GIY-YIG nuclease catalytic core",
      pos_range = c(0.00, 0.65), weight = 3L, evidence_kind = "sequence"
    ),
    "PLD-HKD" = list(
      pattern = "H[A-Z]K.{3,7}D", expected_role = "R", expected_family = "Type II",
      axis = "R", full = "PLD-family HKD motif; two copies and a PLD domain are required",
      pos_range = NULL, weight = 2L, evidence_kind = "sequence"
    ),

    # Type III Res: SF2 motor followed by a C-terminal nuclease domain.
    "ResIII-WA" = list(
      pattern = "(?:[AG].{4}GK[ST]|[GAK]SGK[ST]|AAGK[ST])", expected_role = "R", expected_family = "Type III",
      axis = "R", full = "Type III Res Walker-A/P-loop motif",
      pos_range = c(0.00, 0.75), weight = 2L, evidence_kind = "sequence"
    ),
    "ResIII-WB" = list(
      pattern = "DE.H", expected_role = "R", expected_family = "Type III",
      axis = "R", full = "Type III Res Walker-B DExH motif",
      pos_range = c(0.05, 0.80), weight = 2L, evidence_kind = "sequence"
    ),
    "ResIII-MIII" = list(
      pattern = "[ST]AT", expected_role = "R", expected_family = "Type III",
      axis = "R", full = "Type III Res SF2 motor motif III [ST]AT",
      pos_range = c(0.05, 0.80), weight = 1L, evidence_kind = "sequence"
    ),
    "ResIII-PD" = list(
      pattern = "PD.{1,30}[DE].K", expected_role = "R", expected_family = "Type III",
      axis = "R", full = "Type III Res C-terminal PD-(D/E)xK-like nuclease motif",
      pos_range = c(0.50, 1.00), weight = 3L, evidence_kind = "sequence"
    ),

    # Type IV families are too heterogeneous for one short Mrr regex.
    "TypeIV-profile" = list(
      pattern = NA_character_, expected_role = "R", expected_family = "Type IV",
      axis = "R", full = "Type IV family-specific domain architecture",
      pos_range = NULL, weight = 3L, evidence_kind = "domain_architecture"
    )
  )
}

.dnmb_rebasefinder_motif_pattern <- function(name) {
  definition <- .dnmb_rebasefinder_motif_catalog()[[name]]
  if (base::is.null(definition)) NA_character_ else definition$pattern
}

.dnmb_rebasefinder_motif_all_positions <- function(sequence, pattern) {
  sequence <- .dnmb_rebasefinder_normalize_protein(sequence)
  if (base::is.na(sequence) || base::is.na(pattern) || !base::nzchar(pattern)) {
    return(integer())
  }
  hit <- base::gregexpr(pattern, sequence, perl = TRUE)[[1]]
  if (!base::length(hit) || hit[[1]] < 1L) integer() else base::as.integer(hit)
}

.dnmb_rebasefinder_motif_ordered <- function(position_sets, max_span = Inf) {
  if (!base::length(position_sets) || base::any(base::vapply(position_sets, base::length, integer(1)) == 0L)) {
    return(FALSE)
  }
  walk <- function(i, previous, first) {
    if (i > base::length(position_sets)) return((previous - first) <= max_span)
    candidates <- position_sets[[i]][position_sets[[i]] > previous]
    if (!base::length(candidates)) return(FALSE)
    base::any(base::vapply(candidates, function(candidate) {
      walk(i + 1L, candidate, if (base::is.na(first)) candidate else first)
    }, logical(1)))
  }
  walk(1L, -Inf, NA_real_)
}
