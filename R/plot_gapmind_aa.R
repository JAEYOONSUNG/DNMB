.dnmb_gapmind_aa_nodes <- function() {
  n <- function(id, x, y, label, type = "intermediate") {
    data.frame(id = id, x = x, y = y, label = label, type = type, stringsAsFactors = FALSE)
  }
  s <- 0.5  # uniform step distance
  bx <- 1.0 # backbone x
  # TCA cycle: perfect circle, 8 nodes at 45 deg, center left of backbone
  tca_cx <- 0.0; tca_cy <- 6.5; tca_r <- 1.0  # OAA y = 7.5 = Asp y
  tca_ang <- pi/2 - (0:7) * (2 * pi / 8)  # OAA at 12 o'clock, clockwise
  tx <- tca_cx + tca_r * cos(tca_ang)
  ty <- tca_cy + tca_r * sin(tca_ang)
  # Asp/Glu branch x starts here
  ax <- bx + s  # 1.5
  # Snap TCA y-values to 0.5 grid for downstream pathways
  ay_asp <- round(ty[1] * 2) / 2  # OAA y snapped
  ay_glu <- round(ty[4] * 2) / 2  # a-KG y snapped
  rbind(
    # === Glycolysis backbone (vertical) — G6P stays, G3P+ shift +0.5 ===
    n("G6P",  bx, 13.0, "G6P",      "precursor"),
    n("G3P",  bx, 12.5, "G3P",      "precursor"),
    n("3PG",  bx, 12.0, "3-PG",     "precursor"),
    n("PEP",  bx, 10.5, "PEP",      "precursor"),
    n("PYR",  bx,  9.5, "Pyruvate", "precursor"),
    n("Ala",  bx - s, 9.5, "Ala", "product"),              # Pyruvate → Ala (transamination, left)
    n("AcCoA",bx,  9.0, "AcCoA",    "precursor"),
    # === TCA cycle (perfect circle) ===
    n("OAA",     tx[1], ty[1], "OAA",        "precursor"),
    n("Citrate", tx[2], ty[2], "Citrate",    "precursor"),
    n("Isocit",  tx[3], ty[3], "Isocit", "precursor"),
    n("AKG",     tx[4], ty[4], "a-KG",       "precursor"),
    n("SucCoA",  tx[5], ty[5], "Suc-CoA",    "precursor"),
    n("Succ",    tx[6], ty[6], "Succinate",  "precursor"),
    n("Fum",     tx[7], ty[7], "Fumarate",   "precursor"),
    n("Mal",     tx[8], ty[8], "Malate",     "precursor"),
    # Glyoxylate shunt node (midpoint of Isocit→Mal arc)
    n("Glx",    0.04, 6.6, "Glyoxylate", "precursor"),
    # === PPP branch ===
    n("R5P",  ax, 13.0, "R5P",      "precursor"),
    n("E4P",  ax, 12.5, "E4P",      "precursor"),
    # === Transamination products (snapped to 0.5 grid) ===
    n("Asp",  ax, ay_asp, "Asp", "product"),
    n("Glu",  ax, ay_glu + 0.5, "Glu", "product"),
    # === His: R5P → PRPP → BendHis(horizontal) → His(down) ===
    n("PRPP",    2.0, 13.0, "PRPP", "intermediate"),
    n("BendHis", 5.5, 13.0, "HistolP", "intermediate"),  # hisB product = histidinol-P
    n("His",     5.5, 11.5, "His", "product"),          # hisC,hisN,hisD go down
    # step=s=0.5. Branches: 1-grid diagonal(dx=s,dy=s) or horizontal when possible.
    # === Ser: 3PG +3s→Ser. Gly ↗1s. OAS/Cys horizontal ===
    n("Ser",  bx + 3*s, 12.0, "Ser", "product"),
    n("Gly",  bx + 3*s + s, 12.0 + s, "Gly", "product"),
    n("OAS",  bx + 3*s + s, 12.0, "OAS", "intermediate"),       # horizontal from Ser
    n("Cys",  bx + 3*s + 2*s, 12.0, "Cys", "product"),          # horizontal from OAS
    # === Aromatic: PEP +1s→DAHP +3s→Chor. Prep ↗1s ===
    n("DAHP", bx + 1*s, 10.5, "DAHP",  "intermediate"),
    n("Chor", bx + 4*s, 10.5, "Chor",  "intermediate"),
    # Phe/Tyr: ↗1 diagonal from Chor via Prep
    n("Prep",   bx + 4*s + s, 10.5 + s, "Prep","intermediate"),    # ↗1s diagonal from Chor
    n("Phe",    bx + 4*s + 3*s, 10.5 + s, "Phe", "product"),      # horizontal from Prep
    n("PreTyr", bx + 4*s + 2*s, 10.5 + 2*s, "4HPP","intermediate"),  # 4-hydroxyphenylpyruvate
    n("Tyr",    bx + 4*s + 3*s, 10.5 + 2*s, "Tyr", "product"),   # horizontal from PreTyr, above Phe
    # Trp: horizontal from Chor (no diagonal needed)
    n("Anth", bx + 5*s, 10.5, "Anth","intermediate"),
    n("InGP", bx + 7*s, 10.5, "InGP","intermediate"),
    n("Trp",  bx + 8*s, 10.5, "Trp", "product"),
    # === Pyruvate: PYR +3s→KIV. Val ↗1s. Leu horizontal. Ile ↘1s ===
    n("KIV",  bx + 3*s, 9.5, "KIV",  "intermediate"),
    n("Val",  bx + 3*s + s, 9.5 + s, "Val", "product"),          # ↗1s
    n("KIC",  bx + 3*s + 3*s, 9.5, "KIC", "intermediate"),
    n("Leu",  bx + 3*s + 4*s, 9.5, "Leu", "product"),
    n("2KB",  2.5, 9.0, "2-KB","intermediate"),                # 1x1 diagonal from Thr(3.0,8.5)
    n("Ile",  4.5, 9.0, "Ile", "product"),                    # 4 enzymes × 1 grid each from 2KB
    # === Aspartate: OAA→Asp horizontal. ASA, HSer horizontal. ===
    n("ASA",  ax + 1*s, ay_asp, "ASA",  "intermediate"),
    n("HSer", ax + 2*s, ay_asp, "HSer", "intermediate"),
    n("Thr",  ax + 2*s + s, ay_asp + 2*s, "Thr", "product"),     # ↗ from HSer, 1 grid higher
    n("Met",  ax + 2*s + 5*s, ay_asp, "Met",  "product"),
    n("Asn",  ax + s, ay_asp + 2*s, "Asn", "product"),           # ↗ from Asp, same y as Thr
    n("DHDP", ax + 1*s + s, ay_asp - s, "DHDP","intermediate"),  # ↘1s from ASA
    n("Lys",  ax + 1*s + s + 6*s, ay_asp - s, "Lys","product"),  # horizontal from DHDP
    # === Glutamate: Glu→Gln horizontal. Pro ↘1s. Arg ↘2s ===
    n("Gln",  ax + 1.5, ay_glu + 0.5, "Gln", "product"),         # 0.5 above Glu
    # Pro: place the ProB product (G5P) one row below the arginine branch,
    # then continue horizontally through GSA to Pro
    n("G5P",  ax + 0.5, ay_glu - 0.5, "G5P","intermediate"),
    n("GSA",  ax + 1.0, ay_glu - 0.5, "GSA","intermediate"),
    n("Pro",  ax + 1.5, ay_glu - 0.5, "Pro","product"),
    # Arg: Glu → NAcOrn, then horizontal
    n("NAcOrn",ax + 1.0, ay_glu, "NAcOrn","intermediate"),
    n("Orn",  ax + 1.5, ay_glu, "Orn","intermediate"),
    n("Cit",  ax + 2.0, ay_glu, "Citrulline","intermediate"),
    n("Arg",  ax + 3.0, ay_glu, "Arg","product")
  )
}

# Select best variant for a pathway based on step_status overlap
.dnmb_gapmind_select_variant <- function(pathway_def, pathway_id, step_status) {
  if (is.null(pathway_def$variants)) return(pathway_def)
  ss_steps <- step_status$step_id[step_status$pathway_id == pathway_id & step_status$confidence != "none"]
  best_variant <- NULL
  best_overlap <- -1L
  for (v in pathway_def$variants) {
    all_enzymes <- unlist(v$enzymes)
    overlap <- sum(all_enzymes %in% ss_steps)
    if (overlap > best_overlap) {
      best_overlap <- overlap
      best_variant <- v
    }
  }
  if (is.null(best_variant)) best_variant <- pathway_def$variants[[1]]
  best_variant$family <- pathway_def$family
  best_variant
}

# Pathway route definitions: waypoint nodes + grouped enzyme steps
# Each element of `enzymes` = one reaction segment (may have multiple subunit genes)
# Pathways with `variants` are resolved by .dnmb_gapmind_select_variant()
.dnmb_gapmind_aa_routes <- function(step_status = NULL) {
  raw <- list(
    # breaks = enzymes per waypoint-segment, so waypoint nodes are shared exactly
    his = list(nodes = c("R5P","PRPP","BendHis","His"),
      enzymes = list("prs","hisG","hisE","hisI","hisA","hisF","hisH","hisB","hisC","hisN","hisD"),
      breaks = c(1L, 7L, 3L), family = "Histidine",
      dot_labels = c("PRATP","PRAMP","ProFAR","PRFAR","IGP","Hol-P","Histol","Histal")),
    ser = list(nodes = c("3PG","Ser"),
      enzymes = list("serA","serC","serB"),
      breaks = c(3L), family = "Serine",
      dot_labels = c("3PHP","3PS")),
    gly = list(nodes = c("Ser","Gly"),
      enzymes = list("glyA"),
      breaks = c(1L), family = "Serine"),
    cys = list(nodes = c("Ser","OAS","Cys"),
      enzymes = list("cysE","cysK"),
      breaks = c(1L, 1L), family = "Serine"),
    chorismate = list(nodes = c("PEP","DAHP","Chor"),
      enzymes = list("aroG", c("aroB","aroD"), c("aroE","aroL"), c("aroA","aroC")),
      breaks = c(1L, 3L), family = "Aromatic",
      dot_labels = c("DHS","SHK")),
    phe = list(nodes = c("Chor","Prep","Phe"),
      enzymes = list("cmutase","preph-dehydratase","PPYAT"),
      breaks = c(1L, 2L), family = "Aromatic",
      dot_labels = c("PPA")),
    tyr = list(nodes = c("Chor","Prep","PreTyr","Tyr"),
      enzymes = list("cmutase","pre-dehydr","tyrB"),
      breaks = c(1L, 1L, 1L), family = "Aromatic"),
    trp = list(nodes = c("Chor","Anth","InGP","Trp"),
      enzymes = list("trpE", c("trpD_1","trpD_2"), c("PRAI","IGPS"), c("trpA","trpB")),
      breaks = c(1L, 2L, 1L), family = "Aromatic",
      dot_labels = c("PRAnthr")),
    val = list(nodes = c("PYR","KIV","Val"),
      enzymes = list(c("ilvH","ilvI"), "ilvC", "ilvD", "ilvE"),
      breaks = c(3L, 1L), family = "Pyruvate",
      dot_labels = c("ALAC","DHIV")),
    leu = list(nodes = c("KIV","KIC","Leu"),
      enzymes = list("leuA", c("leuC","leuD"), "leuB", "ilvE"),
      breaks = c(3L, 1L), family = "Pyruvate",
      dot_labels = c("aIPM","bIPM")),
    ile = list(nodes = c("Thr","2KB","Ile"),
      enzymes = list("ilvA", c("ilvH","ilvI"), "ilvC", "ilvD", "ilvE"),
      breaks = c(1L, 4L), family = "Pyruvate",
      dot_labels = c("AHAB","DHMV","KMV")),
    asn = list(nodes = c("Asp","Asn"),
      enzymes = list("aspS2", c("gatC","gatA","gatB")),
      breaks = c(2L), family = "Aspartate",
      dot_labels = c("Asp~tRNA")),
    thr = list(nodes = c("OAA","Asp","ASA","HSer","Thr"),
      enzymes = list("asp-kinase","asd","hom","thrB","thrC"),
      breaks = c(1L, 1L, 1L, 2L), family = "Aspartate",
      dot_labels = c("OPHSer")),
    met = list(
      variants = list(
        transsulfur = list(nodes = c("HSer","Met"),
          enzymes = list("metX","metB","metC","metH","B12-reactivation-domain"),
          breaks = c(5L), dot_labels = c("OAHSer","Cysth","HCys","Met*")),
        direct = list(nodes = c("HSer","Met"),
          enzymes = list("metX","metY","hom","metE"),
          breaks = c(4L), dot_labels = c("OAHSer","HCys","Met*"))
      ), family = "Aspartate"),
    lys = list(
      variants = list(
        dap_dehydr = list(nodes = c("ASA","DHDP","Lys"),
          enzymes = list("dapA","dapB","dapH","dapL","dapX","dapF","lysA"),
          breaks = c(1L, 6L), dot_labels = c("THDP","SDAP","SDAPA","DAP","mDAP")),
        dap_succinyl = list(nodes = c("ASA","DHDP","Lys"),
          enzymes = list("dapA","dapB","dapD","dapC","dapE","dapF","lysA"),
          breaks = c(1L, 6L), dot_labels = c("THDP","SDAP","SDAPA","DAP","mDAP"))
      ), family = "Aspartate"),
    gln = list(nodes = c("Glu","Gln"),
      enzymes = list("gltX", c("gatC","gatA","gatB"), "glnA"),
      breaks = c(3L), family = "Glutamate",
      dot_labels = c("Glu~tRNA","Gln~tRNA")),
    pro = list(
      variants = list(
        ornithine = list(nodes = c("Glu","NAcOrn","GSA","Pro"),
          enzymes = list(c("argJ","argB"), c("argC","argD"), "OAT", "proC"),
          breaks = c(2L, 1L, 1L), dot_labels = c("NAcGlu")),
        direct = list(nodes = c("Glu","G5P","GSA","Pro"),
          enzymes = list("proB","proA","proC"),
          breaks = c(1L, 1L, 1L))
      ), family = "Glutamate"),
    arg = list(nodes = c("Glu","NAcOrn","Orn","Cit","Arg"),
      enzymes = list(c("argJ","argB"), c("argC","argD"), "argI", c("carA","carB"), "argG", "argH"),
      breaks = c(2L, 1L, 1L, 2L), family = "Glutamate",
      dot_labels = c("NAcGlu","ArgSucc"))
  )
  # Resolve variants based on step_status
  resolved <- list()
  for (pid in names(raw)) {
    if (!is.null(raw[[pid]]$variants) && !is.null(step_status)) {
      resolved[[pid]] <- .dnmb_gapmind_select_variant(raw[[pid]], pid, step_status)
    } else if (!is.null(raw[[pid]]$variants)) {
      # No step_status: use first variant
      v1 <- raw[[pid]]$variants[[1]]
      v1$family <- raw[[pid]]$family
      resolved[[pid]] <- v1
    } else {
      resolved[[pid]] <- raw[[pid]]
    }
  }
  resolved
}

# Per-product gene completion stats (found / total enzymes)
.dnmb_gapmind_aa_product_stats <- function(routes, step_status) {
  if (is.null(step_status) || nrow(step_status) == 0) return(NULL)
  ss_key <- paste0(step_status$pathway_id, "::", step_status$step_id)
  stats <- list()
  for (pid in names(routes)) {
    rt <- routes[[pid]]
    enzymes <- rt$enzymes
    if (is.null(enzymes)) enzymes <- as.list(rt$steps)
    all_steps <- unlist(enzymes)
    n_total <- length(all_steps)
    n_found <- sum(vapply(all_steps, function(enz) {
      idx <- match(paste0(pid, "::", enz), ss_key)
      if (is.na(idx)) return(FALSE)
      step_status$confidence[idx] != "none"
    }, logical(1)))
    product_id <- rt$nodes[length(rt$nodes)]
    stats[[pid]] <- data.frame(
      product_id = product_id, pathway_id = pid,
      n_found = n_found, n_total = n_total,
      fraction = if (n_total > 0) n_found / n_total else 0,
      stringsAsFactors = FALSE)
  }
  # Asp/Glu/Ala: trivial transamination (1 step, always present)
  stats[["asp_tx"]] <- data.frame(
    product_id = "Asp", pathway_id = "asp_tx",
    n_found = 1L, n_total = 1L, fraction = 1,
    stringsAsFactors = FALSE)
  stats[["glu_tx"]] <- data.frame(
    product_id = "Glu", pathway_id = "glu_tx",
    n_found = 1L, n_total = 1L, fraction = 1,
    stringsAsFactors = FALSE)
  stats[["ala_tx"]] <- data.frame(
    product_id = "Ala", pathway_id = "ala_tx",
    n_found = 1L, n_total = 1L, fraction = 1,
    stringsAsFactors = FALSE)
  do.call(rbind, stats)
}

# Circle polygon in data coordinates
.dnmb_circle_polygon <- function(cx, cy, r, n_points = 60) {
  theta <- seq(0, 2 * pi, length.out = n_points)
  data.frame(x = cx + r * cos(theta), y = cy + r * sin(theta))
}

# Pie slice polygon (clockwise from 12 o'clock)
.dnmb_pie_polygon <- function(cx, cy, r, fraction, n_points = 60) {
  if (fraction <= 0) return(data.frame(x = numeric(0), y = numeric(0)))
  if (fraction >= 1) return(.dnmb_circle_polygon(cx, cy, r, n_points))
  theta_start <- pi / 2
  theta_end   <- pi / 2 - fraction * 2 * pi
  theta <- seq(theta_start, theta_end, length.out = max(3, round(n_points * fraction)))
  data.frame(x = c(cx, cx + r * cos(theta), cx),
             y = c(cy, cy + r * sin(theta), cy))
}

# ALL metabolite-to-metabolite connections in the network
.dnmb_gapmind_aa_connections <- function(nodes, routes) {
  node_x <- stats::setNames(nodes$x, nodes$id)
  node_y <- stats::setNames(nodes$y, nodes$id)
  # Collect unique directed pairs
  pairs <- list()
  add <- function(a, b) {
    key <- paste0(a, ">>", b)
    pairs[[key]] <<- c(a, b)
  }
  # -- Glycolysis backbone --
  add("G6P","G3P"); add("G3P","3PG")
  add("3PG","PEP"); add("PEP","PYR")
  # -- TCA cycle (circular) --
  add("PYR","AcCoA")
  add("AcCoA","Citrate"); add("OAA","Citrate")
  add("Citrate","Isocit"); add("Isocit","AKG")
  add("AKG","SucCoA"); add("SucCoA","Succ")
  add("Succ","Fum"); add("Fum","Mal"); add("Mal","OAA")
  # -- PPP branch --
  add("G6P","R5P"); add("G3P","E4P")
  # -- Transaminations --
  add("OAA","Asp"); add("AKG","Glu")
  # -- E4P + PEP → DAHP --
  add("E4P","DAHP"); add("PEP","DAHP")
  # -- All route connections --
  for (rt in routes) {
    nds <- rt$nodes
    for (k in seq_along(nds)[-1]) add(nds[k-1], nds[k])
  }
  # -- Cross-pathway connections --
  add("Thr","2KB")   # Ile starts from Thr
  add("Asp","ASA")   # shared by thr/met/lys
  add("ASA","HSer")  # shared by thr/met
  add("ASA","DHDP")  # lys branch
  add("KIV","KIC")   # leu branch from val intermediate
  # Deduplicate (treat A>>B and B>>A as same)
  uniq <- list()
  for (p in pairs) {
    sym <- paste0(min(p[1],p[2]), "::", max(p[1],p[2]))
    uniq[[sym]] <- p
  }
  do.call(rbind, lapply(uniq, function(p) {
    data.frame(from_x = node_x[p[1]], from_y = node_y[p[1]],
               to_x = node_x[p[2]], to_y = node_y[p[2]],
               stringsAsFactors = FALSE)
  }))
}

# Get point at fraction f along polyline
.dnmb_point_on_polyline <- function(wx, wy, f) {
  seg_dx <- diff(wx); seg_dy <- diff(wy)
  seg_len <- sqrt(seg_dx^2 + seg_dy^2)
  cum_len <- c(0, cumsum(seg_len))
  total <- cum_len[length(cum_len)]
  target <- f * total
  si <- max(which(cum_len <= target))
  if (si >= length(cum_len)) si <- length(cum_len) - 1
  local_f <- if (seg_len[si] > 0) (target - cum_len[si]) / seg_len[si] else 0
  c(wx[si] + local_f * seg_dx[si], wy[si] + local_f * seg_dy[si])
}

# Direct line between two points (no elbow, pure diagonal for branches).
# Each enzyme step on a diagonal gets exactly 1 grid cell (dx=0.5, dy=0.5).
.dnmb_elbow_path <- function(x0, y0, x1, y1, arc_n = 8) {
  list(x = c(x0, x1), y = c(y0, y1))
}

# Build enzyme segments using breaks to share waypoint nodes exactly.
# breaks[k] = how many enzyme groups belong to the k-th waypoint segment.
# Waypoint nodes are shared across routes; only intermediate dots are route-specific.
# Returns: list(segments, dots)
.dnmb_build_enzyme_segments <- function(routes, node_x, node_y, step_status) {
  segs <- list(); dots <- list()
  ss_key <- paste0(step_status$pathway_id, "::", step_status$step_id)
  conf_rank <- c(none = 0, low = 1, medium = 2, high = 3)
  for (pid in names(routes)) {
    rt <- routes[[pid]]
    enzymes <- rt$enzymes
    if (is.null(enzymes)) enzymes <- as.list(rt$steps)
    N <- length(enzymes)
    if (N == 0) next
    nds <- rt$nodes
    brks <- rt$breaks
    if (is.null(brks)) brks <- rep(1L, length(nds) - 1)  # fallback: 1 per seg
    enz_idx <- 0
    dot_idx <- 0
    dot_labels_vec <- rt$dot_labels
    for (seg_k in seq_along(brks)) {
      m <- brks[seg_k]  # enzymes in this waypoint segment
      w0 <- nds[seg_k]; w1 <- nds[seg_k + 1]
      x0 <- node_x[w0]; y0 <- node_y[w0]; x1 <- node_x[w1]; y1 <- node_y[w1]
      # Build elbow path
      elbow <- .dnmb_elbow_path(x0, y0, x1, y1)
      ewx <- elbow$x; ewy <- elbow$y
      # m+1 positions along elbow: fraction 0 = start waypoint, 1 = end waypoint
      fracs <- seq(0, 1, length.out = m + 1)
      positions <- lapply(fracs, function(f) {
        p <- .dnmb_point_on_polyline(ewx, ewy, f)
        c(round(p[1] * 2) / 2, round(p[2] * 2) / 2)  # snap to 0.5 grid
      })
      # Use exact node coords for endpoints (don't snap waypoints off-grid)
      positions[[1]] <- c(x0, y0)
      positions[[length(positions)]] <- c(x1, y1)
      for (i in seq_len(m)) {
        enz_idx <- enz_idx + 1
        p0 <- positions[[i]]; p1 <- positions[[i + 1]]
        pm <- c((p0[1]+p1[1])/2, (p0[2]+p1[2])/2)  # midpoint (may be off-grid for diagonals, OK for label)
        enz_group <- if (is.list(enzymes)) enzymes[[enz_idx]] else enzymes[enz_idx]
        confs <- vapply(enz_group, function(enz) {
          m2 <- match(paste0(pid, "::", enz), ss_key)
          if (is.na(m2)) "none" else step_status$confidence[m2]
        }, character(1))
        lts <- vapply(enz_group, function(enz) {
          m2 <- match(paste0(pid, "::", enz), ss_key)
          if (is.na(m2)) NA_character_ else step_status$locus_tag[m2]
        }, character(1))
        best_conf <- names(which.max(conf_rank[confs]))
        lt_found <- lts[!is.na(lts)]
        # Only the FIRST segment of a branching waypoint-pair is curved
        is_branch <- abs(y1 - y0) > 0.01 && abs(x1 - x0) > 0.01 && i == 1
        # Curvature sign: positive = curve up-right, negative = curve down-right
        branch_curv <- if (is_branch) { if (y1 > y0) -0.35 else 0.35 } else 0
        segs[[length(segs) + 1]] <- data.frame(
          x = p0[1], y = p0[2], xend = p1[1], yend = p1[2],
          mx = pm[1], my = pm[2],
          step_label = paste(enz_group, collapse = "\n"),
          pathway_id = pid, family = rt$family,
          confidence = best_conf,
          locus_tag = if (length(lt_found)) paste(lt_found, collapse = "\n") else NA_character_,
          is_branch = is_branch, curvature = branch_curv,
          stringsAsFactors = FALSE)
      }
      # Intermediate dots only (positions 2..m, not start/end waypoints)
      if (m > 1) {
        for (j in 2:m) {
          dot_idx <- dot_idx + 1
          pos <- positions[[j]]
          lbl <- if (!is.null(dot_labels_vec) && dot_idx <= length(dot_labels_vec))
                   dot_labels_vec[dot_idx] else ""
          dots[[length(dots) + 1]] <- data.frame(
            x = pos[1], y = pos[2], pathway_id = pid, family = rt$family,
            label = lbl, stringsAsFactors = FALSE)
        }
      }
    }
  }
  seg_df <- if (length(segs)) do.call(rbind, segs) else data.frame()
  dot_df <- if (length(dots)) do.call(rbind, dots) else data.frame()

  if (nrow(seg_df)) {
    conf_rank <- c(none = 0L, low = 1L, medium = 2L, high = 3L)
    seg_key <- paste(
      round(seg_df$x, 3), round(seg_df$y, 3),
      round(seg_df$xend, 3), round(seg_df$yend, 3),
      seg_df$step_label,
      sep = "|"
    )
    seg_split <- split(seq_len(nrow(seg_df)), seg_key)
    seg_df <- do.call(rbind, lapply(seg_split, function(idx) {
      rows <- seg_df[idx, , drop = FALSE]
      best_idx <- idx[[which.max(conf_rank[rows$confidence])]]
      best <- seg_df[best_idx, , drop = FALSE]
      locus_tags <- unique(stats::na.omit(rows$locus_tag))
      best$locus_tag <- if (length(locus_tags)) paste(locus_tags, collapse = "\n") else NA_character_
      best
    }))
    rownames(seg_df) <- NULL
  }

  if (nrow(dot_df)) {
    dot_key <- paste(
      round(dot_df$x, 3), round(dot_df$y, 3), dot_df$label,
      sep = "|"
    )
    dot_df <- dot_df[!duplicated(dot_key), , drop = FALSE]
    rownames(dot_df) <- NULL
  }

  list(segments = seg_df, dots = dot_df)
}

# Per-step detection status from genbank_table + aa.sum.steps
.dnmb_gapmind_aa_step_status <- function(genbank_table, output_dir = NULL) {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  path_col <- "GapMindAA_pathway_id"; step_col <- "GapMindAA_step_id"
  conf_col <- "GapMindAA_confidence"; bp_col <- "GapMindAA_on_best_path"
  lt_col <- "locus_tag"
  if (!path_col %in% base::names(tbl)) return(NULL)
  tbl <- tbl[!is.na(tbl[[path_col]]) & base::nzchar(tbl[[path_col]]), , drop = FALSE]
  if (!base::nrow(tbl)) return(NULL)
  has_bp <- bp_col %in% names(tbl)
  bp_tbl <- if (has_bp) tbl[!is.na(tbl[[bp_col]]) & tbl[[bp_col]] == TRUE, , drop = FALSE] else tbl
  has_lt <- lt_col %in% names(bp_tbl)
  found <- bp_tbl |>
    dplyr::group_by(pathway_id = .data[[path_col]], step_id = .data[[step_col]]) |>
    dplyr::summarise(
      confidence = dplyr::case_when(
        any(.data[[conf_col]] == "high") ~ "high",
        any(.data[[conf_col]] == "medium") ~ "medium", TRUE ~ "low"),
      locus_tag = if (has_lt) paste(unique(stats::na.omit(.data[[lt_col]])), collapse = "\n") else NA_character_,
      .groups = "drop") |>
    as.data.frame(stringsAsFactors = FALSE)
  # Supplement from aa.sum.steps
  steps_file <- NULL
  if (!is.null(output_dir)) {
    for (f in c(file.path(output_dir, "dnmb_module_gapmindaa", "aa.sum.steps"),
                file.path(output_dir, "dnmb_module_gapmind_aa", "aa.sum.steps"),
                file.path(output_dir, "aa.sum.steps"))) {
      if (file.exists(f)) { steps_file <- f; break }
    }
  }
  if (!is.null(steps_file)) {
    st <- utils::read.delim(steps_file, stringsAsFactors = FALSE)
    if (all(c("onBestPath","pathway","step","score") %in% names(st))) {
      bp_st <- st[st$onBestPath == 1, , drop = FALSE]
      has_locusId <- "locusId" %in% names(bp_st)
      extra <- bp_st |>
        dplyr::group_by(pathway_id = .data$pathway, step_id = .data$step) |>
        dplyr::summarise(
          confidence = dplyr::case_when(
            any(.data$score >= 2) ~ "high",
            any(.data$score >= 1) ~ "medium", TRUE ~ "low"),
          locus_tag = if (has_locusId) paste(unique(stats::na.omit(.data$locusId)), collapse = "\n") else NA_character_,
          .groups = "drop") |>
        as.data.frame(stringsAsFactors = FALSE)
      ekey <- paste0(extra$pathway_id, "::", extra$step_id)
      fkey <- paste0(found$pathway_id, "::", found$step_id)
      found <- rbind(found, extra[!ekey %in% fkey, , drop = FALSE])
    }
  }
  found
}

.dnmb_plot_gapmind_aa_pathway_map <- function(genbank_table, output_dir,
                                               file_stub = "AA_pathway_overview") {
  step_status <- .dnmb_gapmind_aa_step_status(genbank_table, output_dir = output_dir)
  if (is.null(step_status)) return(NULL)

  nodes  <- .dnmb_gapmind_aa_nodes()
  routes <- .dnmb_gapmind_aa_routes(step_status = step_status)
  node_x <- stats::setNames(nodes$x, nodes$id)
  node_y <- stats::setNames(nodes$y, nodes$id)

  # --- Build enzyme segments + intermediate dots ---
  built <- .dnmb_build_enzyme_segments(routes, node_x, node_y, step_status)
  edf  <- built$segments
  ddf  <- built$dots      # intermediate step-boundary circles

  # Edge colors — uniform linewidth, color only distinguishes found/missing
  line_col <- c(high = "#2CA25F", medium = "#FEC44F", low = "#F03B20", none = "#CCCCCC")
  conf_levels <- c("high", "medium", "low", "none")
  edf$confidence <- factor(edf$confidence, levels = conf_levels)
  edf$color <- line_col[as.character(edf$confidence)]
  edf$lw    <- 0.8  # uniform line width
  edf$label_col <- ifelse(edf$confidence == "none", "#AAAAAA", "#222222")

  # Background network connections
  conn <- .dnmb_gapmind_aa_connections(nodes, routes)

  max_x <- max(edf$xend, nodes$x, na.rm = TRUE) + 0.8

  # Node subsets
  products      <- nodes[nodes$type == "product", , drop = FALSE]
  all_inter     <- nodes[nodes$type == "intermediate", , drop = FALSE]
  labeled_inter <- nodes[nodes$type == "intermediate" & nzchar(nodes$label), , drop = FALSE]
  precursors    <- nodes[nodes$type == "precursor", , drop = FALSE]

  lt_df <- edf[!is.na(edf$locus_tag) & nzchar(edf$locus_tag), , drop = FALSE]

  # --- Semantic color gradient ribbons per pathway ---
  # Each pathway gets a unique color; alpha fades from start (0.25) to end (0.05)
  pw_colors <- c(
    his = "#B8860B", ser = "#2171B5", gly = "#6BAED6", cys = "#4292C6",
    chorismate = "#7B68AE", phe = "#9467BD", tyr = "#8C6BB1", trp = "#6A3D9A",
    val = "#E6550D", leu = "#FD8D3C", ile = "#FDAE6B",
    asn = "#74C476", thr = "#31A354", met = "#006D2C", lys = "#41AB5D",
    gln = "#E7298A", pro = "#CE1256", arg = "#DF65B0",
    asp_tx = "#238B45", glu_tx = "#D4507A", ala_tx = "#FF8C00"
  )
  # --- Product pie chart data ---
  pstats <- .dnmb_gapmind_aa_product_stats(routes, step_status)
  pie_product_ids <- if (!is.null(pstats)) unique(pstats$product_id) else character(0)
  products_pie   <- products[products$id %in% pie_product_ids, , drop = FALSE]
  products_nopie <- products[!products$id %in% pie_product_ids, , drop = FALSE]

  pie_r <- 0.09  # match geom_point(size=5) radius in data coords
  pie_list <- list()
  for (i in seq_len(nrow(products_pie))) {
    nid <- products_pie$id[i]
    cx  <- products_pie$x[i]; cy <- products_pie$y[i]
    ps <- pstats[pstats$product_id == nid, , drop = FALSE]
    if (nrow(ps) > 0 && ps$fraction[1] > 0) {
      pie <- .dnmb_pie_polygon(cx, cy, pie_r, ps$fraction[1])
      pie$group <- nid
      pie$fill_color <- pw_colors[ps$pathway_id[1]]
      pie_list[[length(pie_list) + 1]] <- pie
    }
  }
  pie_df_prod <- if (length(pie_list) > 0) do.call(rbind, pie_list) else NULL

  ratio_labels <- if (!is.null(pstats)) {
    rl <- merge(products_pie[, c("id","x","y","label")],
                pstats[, c("product_id","n_found","n_total")],
                by.x = "id", by.y = "product_id", all.x = TRUE)
    rl$ratio_text <- paste0("(", rl$n_found, "/", rl$n_total, ")")
    rl
  } else NULL

  # Ribbon alpha: uniform along each pathway (no left-right gradient).
  # Perpendicular fade (center→edge) is handled by stacked layers below.
  ribbon_segs <- list()
  for (pid in names(routes)) {
    pw_edf <- edf[edf$pathway_id == pid, , drop = FALSE]
    n_seg <- nrow(pw_edf)
    if (n_seg == 0) next
    for (i in seq_len(n_seg)) {
      alpha_i <- 0.18  # constant along the line
      ribbon_segs[[length(ribbon_segs) + 1]] <- data.frame(
        x = pw_edf$x[i], y = pw_edf$y[i],
        xend = pw_edf$xend[i], yend = pw_edf$yend[i],
        color = pw_colors[pid], alpha = alpha_i,
        is_branch = pw_edf$is_branch[i],
        curvature = pw_edf$curvature[i],
        stringsAsFactors = FALSE
      )
    }
  }
  # Add AKG→Glu transamination ribbon
  ribbon_segs[[length(ribbon_segs) + 1]] <- data.frame(
    x = node_x["AKG"], y = node_y["AKG"],
    xend = node_x["Glu"], yend = node_y["Glu"],
    color = pw_colors["glu_tx"], alpha = 0.18,
    is_branch = FALSE, curvature = 0,
    stringsAsFactors = FALSE)
  # Add PYR→Ala transamination ribbon
  ribbon_segs[[length(ribbon_segs) + 1]] <- data.frame(
    x = node_x["PYR"], y = node_y["PYR"],
    xend = node_x["Ala"], yend = node_y["Ala"],
    color = pw_colors["ala_tx"], alpha = 0.18,
    is_branch = FALSE, curvature = 0,
    stringsAsFactors = FALSE)
  ribbon_df <- do.call(rbind, ribbon_segs)

  # Glycolysis backbone (vertical thick line)
  gly_ids <- c("G6P","G3P","3PG","PEP","PYR","AcCoA")
  gly_df <- data.frame(
    x = node_x[gly_ids[-length(gly_ids)]], y = node_y[gly_ids[-length(gly_ids)]],
    xend = node_x[gly_ids[-1]], yend = node_y[gly_ids[-1]],
    stringsAsFactors = FALSE)
  # TCA cycle — smooth circle arc (many points for round look)
  tca_cx <- node_x["OAA"]; tca_cy <- mean(c(node_y["OAA"], node_y["SucCoA"]))
  tca_r  <- (node_y["OAA"] - node_y["SucCoA"]) / 2 + 0.01
  # Recompute center from actual OAA and SucCoA positions for perfect circle
  tca_cx <- (node_x["OAA"] + node_x["SucCoA"]) / 2
  tca_cy <- (node_y["OAA"] + node_y["SucCoA"]) / 2
  tca_r  <- sqrt((node_x["OAA"] - tca_cx)^2 + (node_y["OAA"] - tca_cy)^2)
  theta_seq <- seq(0, 2 * pi, length.out = 100)
  tca_circle <- data.frame(
    x = as.numeric(tca_cx) + tca_r * cos(theta_seq),
    y = as.numeric(tca_cy) + tca_r * sin(theta_seq),
    stringsAsFactors = FALSE)
  # AcCoA → Citrate entry arrow
  ac_cit <- data.frame(x = node_x["AcCoA"], y = node_y["AcCoA"],
    xend = node_x["Citrate"], yend = node_y["Citrate"], stringsAsFactors = FALSE)

  # All nodes table for uniform circle drawing
  all_nodes <- nodes

  # Build base plot, then add ribbon layers
  # Ribbons: center of line = opaque, edges = transparent (gaussian-like via stacked layers)
  p <- ggplot2::ggplot()
  # Draw 4 concentric layers: wide+faint → narrow+opaque
  lw_layers <- c(8, 5.5, 3.5, 1.5)
  alpha_mult <- c(0.25, 0.45, 0.70, 1.0)
  for (li in seq_along(lw_layers)) {
    for (i in seq_len(nrow(ribbon_df))) {
      row <- ribbon_df[i, , drop = FALSE]
      a <- row$alpha * alpha_mult[li]
      if (row$is_branch) {
        p <- p + ggplot2::geom_curve(
          data = row,
          ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
          color = row$color, linewidth = lw_layers[li], alpha = a, curvature = row$curvature)
      } else {
        p <- p + ggplot2::geom_segment(
          data = row,
          ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
          color = row$color, linewidth = lw_layers[li], alpha = a, lineend = "butt")
      }
    }
  }
  p <- p +
    # --- L0a: Glycolysis backbone (light grey) ---
    ggplot2::geom_segment(
      data = gly_df,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#D5D5D5", linewidth = 1.8) +
    # --- L0a2: PPP connections (light grey) ---
    ggplot2::geom_segment(
      data = data.frame(
        x    = c(node_x["G6P"], node_x["G3P"], node_x["E4P"]),
        y    = c(node_y["G6P"], node_y["G3P"], node_y["E4P"]),
        xend = c(node_x["R5P"], node_x["E4P"], node_x["R5P"]),
        yend = c(node_y["R5P"], node_y["E4P"], node_y["R5P"]),
        stringsAsFactors = FALSE),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#D5D5D5", linewidth = 1.2) +
    # --- L0b: TCA cycle (light grey circle) ---
    ggplot2::geom_path(
      data = tca_circle,
      ggplot2::aes(x = .data$x, y = .data$y),
      color = "#D5D5D5", linewidth = 1.5) +
    # --- L0c: AcCoA → Citrate entry ---
    ggplot2::geom_segment(
      data = ac_cit,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#D5D5D5", linewidth = 1.5) +
    # --- L0d: Glyoxylate shunt — Isocit → Mal arc (Glx on midpoint) ---
    ggplot2::geom_curve(
      data = data.frame(x = node_x["Isocit"], y = node_y["Isocit"],
                         xend = node_x["Mal"], yend = node_y["Mal"]),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#B0B0B0", linewidth = 0.7, linetype = "dashed", curvature = -0.3) +
    # --- L0e: AKG → Glu transamination ---
    ggplot2::geom_segment(
      data = data.frame(x = node_x["AKG"], y = node_y["AKG"],
                         xend = node_x["Glu"], yend = node_y["Glu"],
                         stringsAsFactors = FALSE),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#2CA25F", linewidth = 0.8, lineend = "round") +
    # --- L0f: PYR → Ala transamination ---
    ggplot2::geom_segment(
      data = data.frame(x = node_x["PYR"], y = node_y["PYR"],
                         xend = node_x["Ala"], yend = node_y["Ala"],
                         stringsAsFactors = FALSE),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "#2CA25F", linewidth = 0.8, lineend = "round") +
    # --- L1a: Straight enzyme segments (aes color for legend) ---
    ggplot2::geom_segment(
      data = edf[!edf$is_branch, , drop = FALSE],
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend,
                   color = .data$confidence),
      linewidth = 0.8, lineend = "round", show.legend = TRUE) +
    # --- L1b: Curved enzyme segments (fixed color, legend from L1a) ---
    {
      edf_c <- edf[edf$is_branch, , drop = FALSE]
      if (nrow(edf_c) > 0) {
        curve_layers <- lapply(seq_len(nrow(edf_c)), function(j) {
          r <- edf_c[j, , drop = FALSE]
          ggplot2::geom_curve(
            data = r,
            ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
            color = r$color, linewidth = 0.8, curvature = r$curvature,
            show.legend = FALSE)
        })
        curve_layers
      }
    } +
    ggplot2::scale_color_manual(
      name = "GapMind\nConfidence",
      values = c(high = "#2CA25F", medium = "#FEC44F", low = "#F03B20", none = "#CCCCCC"),
      labels = c(high = "High", medium = "Medium", low = "Low", none = "Not found"),
      drop = FALSE) +
    # --- L2: ALL metabolite circles (uniform size) ---
    # Intermediate dots
    ggplot2::geom_point(
      data = ddf, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = "white", color = "#999999", stroke = 0.5) +
    # Intermediate dot metabolite labels (italic, below dot)
    { ddf_labeled <- ddf[nzchar(ddf$label), , drop = FALSE]
      if (nrow(ddf_labeled) > 0) ggplot2::geom_text(
        data = ddf_labeled,
        ggplot2::aes(x = .data$x, y = .data$y - 0.18, label = .data$label),
        size = 1.3, color = "#555555", fontface = "italic")
    } +
    # All intermediates (circles)
    ggplot2::geom_point(
      data = all_inter, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = "white", color = "#999999", stroke = 0.5) +
    # Labeled intermediates (text — italic, below circle)
    ggplot2::geom_text(
      data = labeled_inter,
      ggplot2::aes(x = .data$x, y = .data$y - 0.18, label = .data$label),
      size = 1.5, color = "#444444", fontface = "italic") +
    # Precursors
    ggplot2::geom_point(
      data = precursors, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = "#F0F0F0", color = "#999999", stroke = 0.5) +
    ggplot2::geom_text(
      data = precursors, ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = 1.8, fontface = "bold", color = "#333333") +
    # --- L3: Enzyme + locus_tag — placement by segment orientation ---
    # Vertical up-branch: label LEFT. Vertical down-branch: label RIGHT.
    # Horizontal: label above (down-slope) or below.
    ggplot2::geom_text(
      data = edf,
      ggplot2::aes(
        x = ifelse(abs(.data$xend - .data$x) < 0.05,
                   ifelse(.data$yend > .data$y, .data$mx - 0.15, .data$mx + 0.15),
                   .data$mx),
        y = ifelse(abs(.data$xend - .data$x) < 0.05,
                   .data$my,
                   ifelse(.data$yend < .data$y, .data$my - 0.12, .data$my + 0.10)),
        label = ifelse(is.na(.data$locus_tag) | !nzchar(.data$locus_tag),
                       .data$step_label,
                       paste0(.data$step_label, "\n", .data$locus_tag)),
        hjust = ifelse(abs(.data$xend - .data$x) < 0.05,
                       ifelse(.data$yend > .data$y, 1, 0),
                       0.5),
        vjust = ifelse(abs(.data$xend - .data$x) < 0.05, 0.5,
                       ifelse(.data$yend < .data$y, 1, 0))),
      size = 1.3, color = edf$label_col, fontface = "bold", lineheight = 0.8) +
    # --- L4: Amino acid products (drawn LAST to stay on top) ---
    # White background
    ggplot2::geom_point(
      data = products, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = "white", color = NA, stroke = 0) +
    # Pie fill (fraction of genes found)
    { if (!is.null(pie_df_prod)) ggplot2::geom_polygon(
        data = pie_df_prod,
        ggplot2::aes(x = .data$x, y = .data$y, group = .data$group, fill = .data$fill_color),
        color = NA, show.legend = FALSE) } +
    ggplot2::scale_fill_identity() +
    # Black border (same size as other metabolites)
    ggplot2::geom_point(
      data = products, ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = NA, color = "#111111", stroke = 0.5) +
    # Product name labels
    ggplot2::geom_label(
      data = products, ggplot2::aes(x = .data$x, y = .data$y + 0.20, label = .data$label),
      size = 2.8, fontface = "bold", color = "white", fill = "#333333",
      label.padding = grid::unit(0.12, "lines"), label.size = 0) +
    # Ratio labels inside product circles (white text)
    { if (!is.null(ratio_labels)) ggplot2::geom_text(
        data = ratio_labels,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$ratio_text),
        size = 1.3, color = "white", fontface = "bold") } +
    # --- L5: Chorismate special node (grey fill + white ratio) ---
    ggplot2::geom_point(
      data = data.frame(x = node_x["Chor"], y = node_y["Chor"]),
      ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, size = 5, fill = "#888888", color = "#111111", stroke = 0.5) +
    {
      chor_stats <- pstats[pstats$product_id == "Chor", , drop = FALSE]
      if (nrow(chor_stats) > 0) {
        chor_label <- paste0("(", chor_stats$n_found[1], "/", chor_stats$n_total[1], ")")
        ggplot2::geom_text(
          data = data.frame(x = node_x["Chor"], y = node_y["Chor"]),
          ggplot2::aes(x = .data$x, y = .data$y), label = chor_label,
          size = 1.3, color = "white", fontface = "bold")
      }
    } +
    ggplot2::geom_label(
      data = data.frame(x = node_x["Chor"], y = node_y["Chor"]),
      ggplot2::aes(x = .data$x, y = .data$y + 0.20), label = "Chor",
      size = 2.0, fontface = "bold", color = "white", fill = "#666666",
      label.padding = grid::unit(0.10, "lines"), label.size = 0) +
    # --- Coord & theme ---
    ggplot2::coord_fixed(ratio = 1, xlim = c(-1.3, max_x + 0.5), ylim = c(4.5, 13.5), clip = "off") +
    ggplot2::labs(
      title = "GapMind Amino Acid Biosynthesis Pathway Map",
      subtitle = "Bold border = amino acid product | (n/N) = detected / total genes") +
    # Abbreviation legend as annotation (below plot, wide lines)
    ggplot2::annotate("text", x = -1.7, y = 4.8, hjust = 0, vjust = 1, size = 1.6,
      color = "#555555", lineheight = 1.15, label = paste0(
        "Abbreviations: Ala, alanine; Arg, arginine; Asn, asparagine; Asp, aspartate; Cys, cysteine; Gln, glutamine; Glu, glutamate; Gly, glycine; His, histidine; Ile, isoleucine; Leu, leucine; Lys, lysine;\n",
        "Met, methionine; Phe, phenylalanine; Pro, proline; Ser, serine; Thr, threonine; Trp, tryptophan; Tyr, tyrosine; Val, valine. G6P, glucose-6-P; G3P, glyceraldehyde-3-P; 3-PG, 3-phosphoglycerate;\n",
        "PEP, phosphoenolpyruvate; AcCoA, acetyl-CoA; R5P, ribose-5-P; E4P, erythrose-4-P; OAA, oxaloacetate; a-KG, alpha-ketoglutarate. PRPP, 5-phosphoribosyl-1-pyrophosphate; PRATP/PRAMP, phosphoribosyl-ATP/AMP;\n",
        "ProFAR/PRFAR, phosphoribulosylformimino-AICAR; IGP, imidazoleglycerol-P; Hol-P, imidazole-acetol-P; HistolP, histidinol-P; Histol, histidinol; Histal, histidinal. 3PHP, 3-phosphohydroxypyruvate;\n",
        "3PS, 3-phosphoserine; OAS, O-acetylserine. DAHP, 3-deoxy-D-arabino-heptulosonate-7-P; DHS, dehydroshikimate; SHK, shikimate; Chor, chorismate; Prep, prephenate; PPA, phenylpyruvate; 4HPP, 4-hydroxyphenylpyruvate; Anth, anthranilate;\n",
        "PRAnthr, phosphoribosylanthranilate; InGP, indole-3-glycerol-P. ALAC, acetolactate; DHIV, dihydroxyisovalerate; KIV, ketoisovalerate; aIPM/bIPM, isopropylmalate; KIC, ketoisocaproate; AHAB, acetohydroxybutyrate;\n",
        "DHMV, dihydroxymethylvalerate; KMV, ketomethylvalerate. ASA, aspartate semialdehyde; HSer, homoserine; OPHSer, O-phosphohomoserine; OAHSer, O-acetylhomoserine; Cysth, cystathionine; HCys, homocysteine;\n",
        "DHDP, dihydrodipicolinate; THDP, tetrahydrodipicolinate; SDAP/SDAPA, N-succinyl-DAP intermediates; DAP, LL-diaminopimelate; mDAP, meso-diaminopimelate. NAcGlu, N-acetylglutamate;\n",
        "NAcOrn, N-acetylornithine; GSA, glutamate-semialdehyde; ArgSucc, argininosuccinate")) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12, hjust = 0,
                                         margin = ggplot2::margin(b = 2)),
      plot.subtitle = ggplot2::element_text(size = 7, hjust = 0, color = "#666666",
                                            margin = ggplot2::margin(b = 3)),
      legend.position = c(0.0, 0.92),
      legend.justification = c(0, 1),
      legend.direction = "vertical",
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7),
      legend.key.width = grid::unit(16, "pt"),
      legend.key.height = grid::unit(6, "pt"),
      plot.margin = ggplot2::margin(5, 5, 5, 5))

  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
  # Match PDF dimensions to coord_fixed data range so grid cells look square
  x_range <- diff(range(c(-1.3, max_x + 0.5)))
  y_range <- diff(range(c(4.5, 13.5)))
  plot_scale <- 1.0  # inches per data unit
  .dnmb_module_plot_save(p, pdf_path,
    width = x_range * plot_scale + 1,
    height = y_range * plot_scale + 1)
  list(pdf = pdf_path)
}

.dnmb_plot_gapmind_counts <- function(genbank_table, output_dir, prefix = "GapMindAA", file_stub = "gapmind_aa_pathway_counts", title = "GapMind AA Pathway Counts") {
  tbl <- base::as.data.frame(genbank_table, stringsAsFactors = FALSE)
  path_col <- paste0(prefix, "_pathway_id")
  conf_col <- paste0(prefix, "_confidence")
  if (!path_col %in% base::names(tbl)) {
    return(NULL)
  }
  tbl <- tbl[!is.na(tbl[[path_col]]) & base::nzchar(tbl[[path_col]]), , drop = FALSE]
  if (!base::nrow(tbl)) {
    return(NULL)
  }
  counts <- tbl |>
    dplyr::count(pathway = .data[[path_col]], confidence = if (conf_col %in% base::names(tbl)) .data[[conf_col]] else NA_character_, name = "n") |>
    dplyr::group_by(.data$pathway) |>
    dplyr::mutate(total = sum(.data$n)) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$total), .data$pathway)
  counts$pathway <- factor(counts$pathway, levels = rev(unique(counts$pathway)))
  p <- ggplot2::ggplot(counts, ggplot2::aes(x = .data$n, y = .data$pathway, fill = .data$confidence)) +
    ggplot2::geom_col() +
    ggplot2::labs(
      title = title,
      x = "Gene count",
      y = "Pathway",
      fill = "Confidence"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
  plot_dir <- .dnmb_module_plot_dir(output_dir)
  pdf_path <- file.path(plot_dir, paste0(file_stub, ".pdf"))
  .dnmb_module_plot_save(p, pdf_path, width = 10, height = max(5, 0.25 * length(unique(counts$pathway)) + 2))
  list(pdf = pdf_path)
}

# ======================================================================
# CAZy + Carbon + Membrane Transport Integrated Pathway Map (v2)
# Grid layout: [Extracellular] | Membrane Bilayer | [Cytoplasm] | [Product]
# ======================================================================

# --- CAZy family -> substrate mapping ---
