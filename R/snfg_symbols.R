.dnmb_snfg_star_df <- function(cx, cy, r = 0.15) {
  n <- 5L
  outer_a <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n + 1L)[-(n + 1L)]
  inner_a <- outer_a + pi / n
  ri <- r * 0.52
  data.frame(
    x = c(rbind(r * cos(outer_a) + cx, ri * cos(inner_a) + cx)),
    y = c(rbind(r * sin(outer_a) + cy, ri * sin(inner_a) + cy))
  )
}

# ---- Split diamond (upper triangle colored / lower triangle white) for SNFG uronic acids ----
# This matches the SNFG symbol sheet: top triangle filled, bottom triangle white.
.dnmb_snfg_split_diamond_layers <- function(cx, cy, r = 0.20, fill_color = "#0072BC",
                                            border_color = "black", border_lw = 0.3) {
  # Diamond vertices: top, right, bottom, left
  top   <- c(cx, cy + r)
  right <- c(cx + r, cy)
  bot   <- c(cx, cy - r)
  left  <- c(cx - r, cy)

  # Upper triangle: top-left-right
  upper_df <- data.frame(
    x = c(top[1], right[1], left[1]),
    y = c(top[2], right[2], left[2])
  )

  # Lower triangle: bottom-right-left
  lower_df <- data.frame(
    x = c(bot[1], right[1], left[1]),
    y = c(bot[2], right[2], left[2])
  )

  outline_df <- data.frame(
    x = c(top[1], right[1], bot[1], left[1]),
    y = c(top[2], right[2], bot[2], left[2])
  )
  list(
    ggplot2::geom_polygon(data = lower_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = "white", color = NA, inherit.aes = FALSE),
    ggplot2::geom_polygon(data = upper_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = fill_color, color = NA, inherit.aes = FALSE),
    ggplot2::geom_polygon(data = outline_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = NA, color = border_color, linewidth = border_lw, inherit.aes = FALSE),
    ggplot2::annotate(
      "segment",
      x = left[1], xend = right[1], y = left[2], yend = right[2],
      color = border_color, linewidth = border_lw
    )
  )
}

# ---- Crossed-square polygon for SNFG hexosamine ----
.dnmb_snfg_crossed_square_df <- function(cx, cy, r = 0.15) {
  s <- r * 0.9
  data.frame(
    x = c(cx - s, cx + s, cx + s, cx - s),
    y = c(cy - s, cy - s, cy + s, cy + s)
  )
}

# ---- Diagonally divided square for SNFG Hexosamine ----
# Upper-right / top half is colored, lower-left / bottom half is white.
# The divider runs from top-left to bottom-right (\).
.dnmb_snfg_crossed_square_layers <- function(cx, cy, s = 0.15, fill_color = "#0072BC",
                                             border_color = "black", border_lw = 0.3) {
  lt <- c(cx - s, cy + s)
  rt <- c(cx + s, cy + s)
  rb <- c(cx + s, cy - s)
  lb <- c(cx - s, cy - s)

  upper_df <- data.frame(
    x = c(lt[1], rt[1], rb[1]),
    y = c(lt[2], rt[2], rb[2])
  )
  lower_df <- data.frame(
    x = c(lt[1], rb[1], lb[1]),
    y = c(lt[2], rb[2], lb[2])
  )
  outline_df <- data.frame(
    x = c(lt[1], rt[1], rb[1], lb[1]),
    y = c(lt[2], rt[2], rb[2], lb[2])
  )

  list(
    ggplot2::geom_polygon(
      data = lower_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = "white", color = NA, inherit.aes = FALSE
    ),
    ggplot2::geom_polygon(
      data = upper_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = fill_color, color = NA, inherit.aes = FALSE
    ),
    ggplot2::geom_polygon(
      data = outline_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = NA, color = border_color, linewidth = border_lw, inherit.aes = FALSE
    ),
    ggplot2::annotate(
      "segment",
      x = lt[1], xend = rb[1], y = lt[2], yend = rb[2],
      color = border_color, linewidth = border_lw
    )
  )
}

# ---- Pentagon polygon for SNFG ketose/assigned (Fructose, Sorbitol, etc.) ----
.dnmb_snfg_pentagon_df <- function(cx, cy, r = 0.18) {
  angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = 6L)[-(6L)]
  data.frame(x = r * cos(angles) + cx, y = r * sin(angles) + cy)
}

# ---- Divided triangle (upper white / lower colored) for SNFG DeoxyhexNAc ----
.dnmb_snfg_divided_triangle_layers <- function(cx, cy, r = 0.20, fill_color = "#0072BC") {
  # Equilateral-ish triangle: top vertex, bottom-left, bottom-right
  top <- c(cx, cy + r)
  bl  <- c(cx - r, cy - r * 0.7)
  br  <- c(cx + r, cy - r * 0.7)
  mid_l <- c((top[1] + bl[1]) / 2, (top[2] + bl[2]) / 2)
  mid_r <- c((top[1] + br[1]) / 2, (top[2] + br[2]) / 2)
  # Upper half: top, mid_r, mid_l → white
  upper_df <- data.frame(x = c(top[1], mid_r[1], mid_l[1]),
                          y = c(top[2], mid_r[2], mid_l[2]))
  # Lower half: mid_l, mid_r, br, bl → colored
  lower_df <- data.frame(x = c(mid_l[1], mid_r[1], br[1], bl[1]),
                          y = c(mid_l[2], mid_r[2], br[2], bl[2]))
  # Full outline
  outline_df <- data.frame(x = c(top[1], br[1], bl[1]),
                            y = c(top[2], br[2], bl[2]))
  list(
    ggplot2::geom_polygon(data = lower_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = fill_color, color = NA, inherit.aes = FALSE),
    ggplot2::geom_polygon(data = upper_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = "white", color = NA, inherit.aes = FALSE),
    ggplot2::geom_polygon(data = outline_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill = NA, color = "black", linewidth = 0.2, inherit.aes = FALSE),
    # Dividing line
    ggplot2::geom_segment(data = data.frame(x = mid_l[1], xend = mid_r[1],
      y = mid_l[2], yend = mid_r[2]),
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
      color = "black", linewidth = 0.2, inherit.aes = FALSE)
  )
}

# ---- Flat diamond (wider than tall) for SNFG 3,9-dideoxy-nonulosonic ----
.dnmb_snfg_flat_diamond_df <- function(cx, cy, r = 0.20) {
  hw <- r * 1.2  # wider
  hh <- r * 0.6  # shorter
  data.frame(
    x = c(cx, cx + hw, cx, cx - hw),
    y = c(cy + hh, cy, cy - hh, cy)
  )
}

# ---- SNFG official color palette ----
.dnmb_snfg_colors <- function() {
  list(
    white      = "#FFFFFF",
    blue       = "#0072BC",
    green      = "#05A551",
    yellow     = "#FFD400",
    orange     = "#F47920",
    pink       = "#F59EA1",
    purple     = "#A54399",
    light_blue = "#8FCCE9",
    brown      = "#A17A4D",
    red        = "#ED1C25"
  )
}

# ---- Border color helper: darker shade of the SNFG fill ----
.dnmb_snfg_border_color <- function(fill_color) {
  border_map <- c(
    "#FFFFFF" = "#BDBDBD",
    "#0072BC" = "#004A7C",
    "#05A551" = "#006B34",
    "#FFD400" = "#CCA800",
    "#F47920" = "#C55A10",
    "#F59EA1" = "#C47078",
    "#A54399" = "#7E2F76",
    "#8FCCE9" = "#5E9EBB",
    "#A17A4D" = "#6F5332",
    "#ED1C25" = "#B51219"
  )
  if (fill_color %in% names(border_map)) return(unname(border_map[fill_color]))
  rgb <- grDevices::col2rgb(fill_color) / 255
  darker <- pmax(0, rgb * 0.65)
  grDevices::rgb(darker[1], darker[2], darker[3])
}

# ---- SNFG full lookup table ----
# shape_type: 'circle'=21L, 'square'=22L, 'diamond'=23L, 'triangle'=24L,
#   custom: 'star'=-1L, 'crossed_square'=-2L, 'pentagon'=-3L,
#           'divided_diamond'=23L (rendered as split), 'flat_rectangle'=-4L,
#           'hexagon'=-5L, 'filled_diamond'=23L, 'divided_triangle'=-6L
.dnmb_snfg_table <- function() {
  cl <- .dnmb_snfg_colors()
  # 10-color order per row: white, blue, green, yellow, orange, pink, purple, lightblue, brown, red
  color_cycle <- c(cl$white, cl$blue, cl$green, cl$yellow, cl$orange,
                   cl$pink, cl$purple, cl$light_blue, cl$brown, cl$red)
  data.frame(
    canonical = c(
      # Hexose (Filled Circle) — 8 of 10 columns used
      "Glc","Man","Gal","Gul","Alt","All","Tal","Ido",
      # HexNAc (Filled Square)
      "GlcNAc","ManNAc","GalNAc","GulNAc","AltNAc","AllNAc","TalNAc","IdoNAc",
      # Hexosamine (Crossed Square)
      "GlcN","ManN","GalN","GulN","AltN","AllN","TalN","IdoN",
      # Hexuronate (Divided Diamond)
      "GlcA","ManA","GalA","GulA","AltA","AllA","TalA","IdoA",
      # Deoxyhexose (Filled Triangle)
      "Qui","Rha","6dGul","6dAlt","6dTal","Fuc",
      # DeoxyhexNAc (Divided Triangle)
      "QuiNAc","RhaNAc","6dAltNAc","6dTalNAc","FucNAc",
      # Di-deoxyhexose (Flat Rectangle)
      "Oli","Tyv","Abe","Par","Dig","Col",
      # Pentose (Filled Star)
      "Ara","Lyx","Xyl","Rib",
      # Nonulosonate (Filled Diamond)
      "Kdn","Neu5Ac","Neu5Gc","Neu","Sia",
      # 3,9-dideoxy-nonulosonic (Flat Diamond)
      "Pse","Leg","Aci","4eLeg",
      # Unknown/Other (Flat Hexagon)
      "Bac","LDmanHep","Kdo","Dha","DDmanHep","MurNAc","MurNGc","Mur",
      # Assigned/Ketose (Pentagon)
      "Api","Fru","Tag","Sor","Psi"
    ),
    shape_code = c(
      rep(21L, 8),   # circle (hexose)
      rep(22L, 8),   # square (HexNAc)
      rep(-2L, 8),   # crossed square (hexosamine)
      rep(23L, 8),   # divided diamond (hexuronate) — rendered as split
      rep(24L, 6),   # triangle (deoxyhexose)
      rep(-6L, 5),   # divided triangle (DeoxyhexNAc)
      rep(-4L, 6),   # flat rectangle (di-deoxyhexose)
      rep(-1L, 4),   # star (pentose)
      rep(23L, 5),   # filled diamond (nonulosonate) — NOT split
      rep(-7L, 4),   # flat diamond (3,9-dideoxy-nonulosonic)
      rep(-5L, 8),   # hexagon (unknown/bacterial)
      rep(-3L, 5)    # pentagon (ketose)
    ),
    is_split = c(
      rep(FALSE, 8), rep(FALSE, 8), rep(FALSE, 8),
      rep(TRUE, 8),   # hexuronate = split diamond
      rep(FALSE, 6),
      rep(TRUE, 5),   # DeoxyhexNAc = split triangle
      rep(FALSE, 6), rep(FALSE, 4),
      rep(FALSE, 5),  # nonulosonate = filled diamond (not split)
      rep(FALSE, 4),  # flat diamond
      rep(FALSE, 8),
      rep(FALSE, 5)
    ),
    # Column assignments per NCBI table:
    # col: 1=white, 2=blue, 3=green, 4=yellow, 5=orange, 6=pink, 7=purple, 8=lightblue, 9=brown, 10=red
    fill_color = c(
      # Hexose: Glc=blue(2), Man=green(3), Gal=yellow(4), Gul=orange(5), Alt=pink(6), All=purple(7), Tal=lightblue(8), Ido=brown(9)
      cl$blue, cl$green, cl$yellow, cl$orange, cl$pink, cl$purple, cl$light_blue, cl$brown,
      # HexNAc: same color order
      cl$blue, cl$green, cl$yellow, cl$orange, cl$pink, cl$purple, cl$light_blue, cl$brown,
      # Hexosamine: same color order
      cl$blue, cl$green, cl$yellow, cl$orange, cl$pink, cl$purple, cl$light_blue, cl$brown,
      # Hexuronate: same color order
      cl$blue, cl$green, cl$yellow, cl$orange, cl$pink, cl$purple, cl$light_blue, cl$brown,
      # Deoxyhexose: Qui=blue(2), Rha=green(3), 6dGul=orange(5), 6dAlt=pink(6), 6dTal=lightblue(8), Fuc=red(10)
      cl$blue, cl$green, cl$orange, cl$pink, cl$light_blue, cl$red,
      # DeoxyhexNAc: QuiNAc=blue(2), RhaNAc=green(3), 6dAltNAc=yellow(4), 6dTalNAc=orange(5), FucNAc=red(10)
      cl$blue, cl$green, cl$yellow, cl$orange, cl$red,
      # Di-deoxyhexose: Oli=blue(2), Tyv=green(3), Abe=orange(5), Par=pink(6), Dig=purple(7), Col=lightblue(8)
      cl$blue, cl$green, cl$orange, cl$pink, cl$purple, cl$light_blue,
      # Pentose: Ara=blue(2), Lyx=green(3), Xyl=yellow(4), Rib=orange(5)
      cl$blue, cl$green, cl$yellow, cl$orange,
      # Nonulosonate: Kdn=green(3), Neu5Ac=purple(7), Neu5Gc=lightblue(8), Neu=brown(9), Sia=red(10)
      cl$green, cl$purple, cl$light_blue, cl$brown, cl$red,
      # Flat Diamond: Pse=green(3), Leg=yellow(4), Aci=orange(5), 4eLeg=lightblue(8)
      cl$green, cl$yellow, cl$orange, cl$light_blue,
      # Hexagon: Bac=blue(2), LDmanHep=green(3), Kdo=yellow(4), Dha=orange(5), DDmanHep=pink(6), MurNAc=purple(7), MurNGc=lightblue(8), Mur=brown(9)
      cl$blue, cl$green, cl$yellow, cl$orange, cl$pink, cl$purple, cl$light_blue, cl$brown,
      # Pentagon: Api=blue(2), Fru=green(3), Tag=yellow(4), Sor=orange(5), Psi=pink(6)
      cl$blue, cl$green, cl$yellow, cl$orange, cl$pink
    ),
    # NCBI column index (1-10) for grid alignment
    color_col = c(
      # Hexose: cols 2-9
      2L,3L,4L,5L,6L,7L,8L,9L,
      # HexNAc: cols 2-9
      2L,3L,4L,5L,6L,7L,8L,9L,
      # Hexosamine: cols 2-9
      2L,3L,4L,5L,6L,7L,8L,9L,
      # Hexuronate: cols 2-9
      2L,3L,4L,5L,6L,7L,8L,9L,
      # Deoxyhexose: Qui=2, Rha=3, 6dGul=5, 6dAlt=6, 6dTal=8, Fuc=10
      2L,3L,5L,6L,8L,10L,
      # DeoxyhexNAc: QuiNAc=2, RhaNAc=3, 6dAltNAc=4, 6dTalNAc=5, FucNAc=10
      2L,3L,4L,5L,10L,
      # Di-deoxyhexose: Oli=2, Tyv=3, Abe=5, Par=6, Dig=7, Col=8
      2L,3L,5L,6L,7L,8L,
      # Pentose: Ara=2, Lyx=3, Xyl=4, Rib=5
      2L,3L,4L,5L,
      # Nonulosonate: Kdn=3, Neu5Ac=7, Neu5Gc=8, Neu=9, Sia=10
      3L,7L,8L,9L,10L,
      # Flat Diamond: Pse=3, Leg=4, Aci=5, 4eLeg=8
      3L,4L,5L,8L,
      # Hexagon: Bac=2, LDmanHep=3, Kdo=4, Dha=5, DDmanHep=6, MurNAc=7, MurNGc=8, Mur=9
      2L,3L,4L,5L,6L,7L,8L,9L,
      # Pentagon: Api=2, Fru=3, Tag=4, Sor=5, Psi=6
      2L,3L,4L,5L,6L
    ),
    stringsAsFactors = FALSE
  )
}

# ---- SNFG synonym resolver ----
.dnmb_snfg_resolve_synonym <- function(name) {
  synonyms <- list(
    Glc = c("glucose","Glucose","D-glucose","Glc","glc"),
    Gal = c("galactose","Galactose","D-galactose","Gal","gal"),
    Man = c("mannose","Mannose","D-mannose","Man","man"),
    Fru = c("fructose","Fructose","D-fructose","Fru","fru"),
    GlcNAc = c("NAG","GlcNAc","N-acetylglucosamine","N-Acetylglucosamine","glcnac","nag"),
    GalNAc = c("GalNAc","N-acetylgalactosamine","galnac"),
    ManNAc = c("ManNAc","N-acetylmannosamine","mannac"),
    GlcN = c("glucosamine","Glucosamine","GlcN","glcn"),
    GlcA = c("glucuronate","Glucuronate","glucuronic acid","GlcA","glca","gluconate"),
    GalA = c("galacturonate","Galacturonate","galacturonic acid","GalA","gala"),
    Fuc = c("fucose","Fucose","L-fucose","Fuc","fuc"),
    Rha = c("rhamnose","Rhamnose","L-rhamnose","Rha","rha"),
    Xyl = c("xylose","Xylose","D-xylose","Xyl","xyl"),
    Ara = c("arabinose","Arabinose","L-arabinose","Ara","ara"),
    Rib = c("ribose","Ribose","D-ribose","Rib","rib"),
    MurNAc = c("MurNAc","N-acetylmuramic acid","muramic acid"),
    Neu5Ac = c("Neu5Ac","sialic acid","NeuNAc","NANA"),
    Sor = c("sorbitol","Sorbitol","glucitol","Sor"),
    Tag = c("tagatose","Tagatose","Tag"),
    Gul = c("gulose","Gulose","D-gulose","Gul","gul"),
    Alt = c("altrose","Altrose","D-altrose","Alt","alt"),
    All = c("allose","Allose","D-allose","All"),
    Tal = c("talose","Talose","D-talose","Tal","tal"),
    Ido = c("idose","Idose","D-idose","Ido","ido"),
    Qui = c("quinovose","Qui"),
    Kdo = c("Kdo","KDO","3-deoxy-D-manno-octulosonate"),
    Maltose = c("maltose","Maltose"),
    Cellobiose = c("cellobiose","Cellobiose"),
    Trehalose = c("trehalose","Trehalose"),
    Sucrose = c("sucrose","Sucrose"),
    Lactose = c("lactose","Lactose"),
    Raffinose = c("raffinose","Raffinose"),
    Chitobiose = c("chitobiose","Chitobiose"),
    Melibiose = c("melibiose","Melibiose"),
    Amylopectin = c("amylopectin","Amylopectin"),
    Chitosan = c("chitosan","Chitosan"),
    Galactomannan = c("galactomannan","Galactomannan"),
    Rhamnogalacturonan = c("rhamnogalacturonan","Rhamnogalacturonan","RG-I","RG-II"),
    Levan = c("levan","Levan"),
    Dextran = c("dextran","Dextran"),
    Pullulan = c("pullulan","Pullulan"),
    Hyaluronan = c("hyaluronan","Hyaluronan","hyaluronic acid","Hyaluronic acid"),
    Chondroitin = c("chondroitin","Chondroitin","chondroitin sulfate","Chondroitin sulfate"),
    Heparin = c("heparin","Heparin"),
    Isomaltose = c("isomaltose","Isomaltose"),
    IdoA = c("iduronate","Iduronate","iduronic acid","IdoA","idoa"),
    Amylose = c("amylose","Amylose"),
    Curdlan = c("curdlan","Curdlan"),
    Laminarin = c("laminarin","Laminarin","laminaran"),
    Agarose = c("agarose","Agarose"),
    Alginate = c("alginate","Alginate","alginic acid"),
    Glucuronoxylan = c("glucuronoxylan","Glucuronoxylan"),
    Glucomannan = c("glucomannan","Glucomannan","konjac glucomannan"),
    Xyloglucan = c("xyloglucan","Xyloglucan"),
    `Heparan sulfate` = c("heparan sulfate","Heparan sulfate","heparan","Heparan"),
    `Dermatan sulfate` = c("dermatan sulfate","Dermatan sulfate","dermatan","Dermatan"),
    `Keratan sulfate` = c("keratan sulfate","Keratan sulfate","keratan","Keratan"),
    `beta-Glucan (1-3/1-6)` = c("beta-Glucan (1-3/1-6)","beta-glucan","Beta-glucan","beta-Glucan"),
    ManA = c("mannuronate","Mannuronate","mannuronic acid","ManA","mana"),
    GulA = c("guluronate","Guluronate","guluronic acid","GulA","gula"),
    `3,6-AG` = c("3,6-anhydrogalactose","3,6-AG","3,6-AnGal","anhydrogalactose")
  )
  name_lc <- tolower(name)
  for (canonical in names(synonyms)) {
    if (name_lc %in% tolower(synonyms[[canonical]])) return(canonical)
  }
  name
}

# ---- SNFG oligosaccharide / polysaccharide definitions ----
.dnmb_snfg_oligosaccharide <- function() {
  list(
    Maltose    = list(units = c("Glc","Glc"), bond = "alpha-1,4", linetype = "solid"),
    Cellobiose = list(units = c("Glc","Glc"), bond = "beta-1,4",  linetype = "dashed"),
    Trehalose  = list(units = c("Glc","Glc"), bond = "alpha-1,1", linetype = "solid"),
    Sucrose    = list(units = c("Glc","Fru"), bond = "alpha-1,2", linetype = "solid"),
    Lactose    = list(units = c("Gal","Glc"), bond = "beta-1,4",  linetype = "dashed"),
    Raffinose  = list(units = c("Gal","Glc","Fru"), bond = "alpha-1,6/alpha-1,2", linetype = "solid"),
    Chitobiose = list(units = c("GlcNAc","GlcNAc"), bond = "beta-1,4", linetype = "dashed"),
    Melibiose  = list(units = c("Gal","Glc"), bond = "alpha-1,6", linetype = "solid"),
    Starch     = list(units = c("Glc","Glc","Glc","Glc","Glc"), bond = "alpha-1,4", linetype = "solid"),
    Cellulose  = list(units = c("Glc","Glc","Glc","Glc","Glc"), bond = "beta-1,4", linetype = "dashed"),
    Chitin     = list(units = c("GlcNAc","GlcNAc","GlcNAc","GlcNAc"), bond = "beta-1,4", linetype = "dashed"),
    Xylan      = list(units = c("Xyl","Xyl","Xyl","Xyl"), bond = "beta-1,4", linetype = "dashed"),
    Mannan     = list(units = c("Man","Man","Man","Man"), bond = "beta-1,4", linetype = "dashed"),
    Pectin     = list(units = c("GalA","GalA","GalA","GalA"), bond = "alpha-1,4", linetype = "solid"),
    Inulin     = list(units = c("Fru","Fru","Fru","Fru"), bond = "beta-2,1", linetype = "dashed"),
    Peptidoglycan = list(units = c("GlcNAc","MurNAc","GlcNAc","MurNAc"), bond = "beta-1,4", linetype = "dashed"),
    Amylopectin = list(units = c("Glc","Glc","Glc","Glc","Glc"), bond = "alpha-1,4/alpha-1,4/alpha-1,6/alpha-1,4", linetype = "solid"),
    Chitosan    = list(units = c("GlcN","GlcN","GlcN","GlcN"), bond = "beta-1,4", linetype = "dashed"),
    Galactomannan = list(units = c("Man","Man","Man","Man"), bond = "beta-1,4", linetype = "dashed"),
    Rhamnogalacturonan = list(units = c("GalA","Rha","GalA","Rha"), bond = "alpha-1,2/alpha-1,4/alpha-1,2", linetype = "solid"),
    Levan       = list(units = c("Fru","Fru","Fru","Fru"), bond = "beta-2,6", linetype = "dashed"),
    Dextran     = list(units = c("Glc","Glc","Glc","Glc"), bond = "alpha-1,6", linetype = "solid"),
    Pullulan    = list(units = c("Glc","Glc","Glc","Glc","Glc","Glc"), bond = "alpha-1,4/alpha-1,4/alpha-1,6/alpha-1,4/alpha-1,4", linetype = "solid"),
    Hyaluronan  = list(units = c("GlcA","GlcNAc","GlcA","GlcNAc"), bond = "beta-1,3/beta-1,4/beta-1,3", linetype = "dashed"),
    Chondroitin = list(units = c("GlcA","GalNAc","GlcA","GalNAc"), bond = "beta-1,3/beta-1,4/beta-1,3", linetype = "dashed"),
    Heparin     = list(units = c("IdoA","GlcNAc","IdoA","GlcNAc"), bond = "alpha-1,4/alpha-1,4/alpha-1,4", linetype = "solid"),
    Isomaltose  = list(units = c("Glc","Glc"), bond = "alpha-1,6", linetype = "solid"),
    Amylose     = list(units = c("Glc","Glc","Glc","Glc","Glc"), bond = "alpha-1,4", linetype = "solid"),
    Curdlan     = list(units = c("Glc","Glc","Glc","Glc"), bond = "beta-1,3", linetype = "dashed"),
    Laminarin   = list(units = c("Glc","Glc","Glc","Glc"), bond = "beta-1,3/beta-1,3/beta-1,6", linetype = "dashed"),
    Agarose     = list(units = c("Gal","Gal","Gal","Gal"), bond = "beta-1,4/alpha-1,3/beta-1,4", linetype = "dashed"),
    Alginate    = list(units = c("ManA","GulA","ManA","GulA"), bond = "beta-1,4", linetype = "dashed"),
    Glucuronoxylan = list(units = c("Xyl","Xyl","Xyl","GlcA"), bond = "beta-1,4/beta-1,4/alpha-1,2", linetype = "dashed"),
    Glucomannan = list(units = c("Man","Glc","Man","Glc"), bond = "beta-1,4", linetype = "dashed"),
    Xyloglucan  = list(units = c("Glc","Glc","Glc","Xyl"), bond = "beta-1,4/beta-1,4/alpha-1,6", linetype = "dashed"),
    `Heparan sulfate` = list(units = c("GlcA","GlcNAc","GlcA","GlcNAc"), bond = "beta-1,4/alpha-1,4/beta-1,4", linetype = "dashed"),
    `Dermatan sulfate` = list(units = c("IdoA","GalNAc","IdoA","GalNAc"), bond = "alpha-1,3/beta-1,4/alpha-1,3", linetype = "solid"),
    `Keratan sulfate` = list(units = c("Gal","GlcNAc","Gal","GlcNAc"), bond = "beta-1,4/beta-1,3/beta-1,4", linetype = "dashed"),
    `beta-Glucan (1-3/1-6)` = list(units = c("Glc","Glc","Glc","Glc"), bond = "beta-1,3/beta-1,3/beta-1,6", linetype = "dashed")
  )
}

# ---- Unified SNFG lookup (replaces old .dnmb_snfg_style + .dnmb_snfg_lookup) ----
# Accepts raw sugar_type strings (synonyms resolved automatically)
# Returns data.frame with columns: snfg_shape, snfg_color, snfg_fill
.dnmb_snfg_lookup <- function(sugar_types) {
  tbl <- .dnmb_snfg_table()
  out <- data.frame(
    snfg_shape = integer(length(sugar_types)),
    snfg_color = character(length(sugar_types)),
    snfg_fill  = character(length(sugar_types)),
    stringsAsFactors = FALSE)
  for (i in seq_along(sugar_types)) {
    canonical <- .dnmb_snfg_resolve_synonym(sugar_types[i])
    idx <- match(canonical, tbl$canonical)
    if (is.na(idx)) {
      out$snfg_shape[i] <- 21L
      out$snfg_color[i] <- "#AAAAAA"
      out$snfg_fill[i]  <- "#DDDDDD"
    } else {
      out$snfg_shape[i] <- tbl$shape_code[idx]
      out$snfg_color[i] <- tbl$fill_color[idx]
      out$snfg_fill[i]  <- if (tbl$is_split[idx]) "#FFFFFF" else tbl$fill_color[idx]
    }
  }
  out
}

# ---- SNFG metabolite composition (uses oligosaccharide table + synonym resolver) ----
.dnmb_snfg_metabolite_composition <- function(metabolite) {
  canonical <- .dnmb_snfg_resolve_synonym(metabolite)
  oligo <- .dnmb_snfg_oligosaccharide()
  # Check if it's an oligosaccharide

  if (grepl("\\+", canonical)) {
    parts <- trimws(strsplit(canonical, "\\+")[[1]])
    parts <- parts[nzchar(parts)]
    if (length(parts) > 0) {
      monomers <- vapply(parts, .dnmb_snfg_resolve_synonym, character(1))
      return(list(monomers = monomers, linetype = NA, bond_labels = NA_character_,
                  phospho = FALSE, is_sugar = TRUE))
    }
  }

  if (canonical %in% names(oligo)) {
    o <- oligo[[canonical]]
    # Parse bond description into display label(s) and linetype
    bond_parts <- strsplit(o$bond, "/")[[1]]
    bond_labels <- vapply(bond_parts, function(b) {
      b <- trimws(b)
      b <- sub("^alpha", "alpha", b)
      b <- sub("^beta",  "beta", b)
      b
    }, character(1), USE.NAMES = FALSE)
    return(list(monomers = o$units, linetype = o$linetype, bond_labels = bond_labels,
                phospho = FALSE, is_sugar = TRUE))
  }
  # Check if it's a known monosaccharide
  tbl <- .dnmb_snfg_table()
  if (canonical %in% tbl$canonical) {
    return(list(monomers = canonical, linetype = NA, bond_labels = NA_character_,
                phospho = FALSE, is_sugar = TRUE))
  }
  # Non-sugar metabolites
  nonsugars <- c("ethanol","acetate","L-lactate","D-lactate","citrate",
                 "succinate","fumarate","L-malate","glycerol","mannitol","myoinositol")
  if (tolower(metabolite) %in% tolower(nonsugars)) {
    return(list(monomers = character(0), linetype = NA, bond_labels = NA_character_,
                phospho = FALSE, is_sugar = FALSE))
  }
  # Default: try as single monosaccharide by synonym
  list(monomers = canonical, linetype = NA, bond_labels = NA_character_,
       phospho = FALSE, is_sugar = TRUE)
}

# ---- Render a single SNFG symbol at given center (cx, cy) ----
# All shapes unified to HexNAc square as reference:
#   r = half-side of the reference square. The square's bounding circle = r * sqrt(2).
#   All other shapes' bounding circles match this same radius = r * sqrt(2).
.dnmb_snfg_render_symbol <- function(cx, cy, sugar_type, r = 0.33,
                                      point_size = 9, active = TRUE) {
  canonical <- .dnmb_snfg_resolve_synonym(sugar_type)
  tbl <- .dnmb_snfg_table()
  idx <- match(canonical, tbl$canonical)
  blw <- 0.3
  layers <- list()

  if (is.na(idx)) {
    bdr <- .dnmb_snfg_border_color("#AAAAAA")
    layers[[1]] <- ggplot2::annotate("point", x = cx, y = cy,
      shape = 21L, size = point_size, fill = "#AAAAAA", color = bdr, stroke = blw)
    return(layers)
  }

  pch <- tbl$shape_code[idx]
  col <- tbl$fill_color[idx]
  is_split <- tbl$is_split[idx]
  bdr <- .dnmb_snfg_border_color(col)

  # s = square half-side (reference). R = bounding circle radius = s * sqrt(2)
  s <- r
  R <- s * sqrt(2)

  # Helper: polygon layer
  gpoly <- function(df, fill, col = bdr, lw = blw)
    ggplot2::geom_polygon(data = df, ggplot2::aes(x = .data$x, y = .data$y),
      fill = fill, color = col, linewidth = lw, inherit.aes = FALSE)
  gpoly_na <- function(df, fill)
    ggplot2::geom_polygon(data = df, ggplot2::aes(x = .data$x, y = .data$y),
      fill = fill, color = NA, inherit.aes = FALSE)

  if (pch == 21L) {
    # Circle (Hexose): radius = R (same bounding circle as square)
    circ_angles <- seq(0, 2 * pi, length.out = 61)[-61]
    circ_df <- data.frame(x = R * cos(circ_angles) + cx, y = R * sin(circ_angles) + cy)
    layers[[1]] <- gpoly(circ_df, col)

  } else if (pch == 22L) {
    # Square (HexNAc): half-side = s
    sq_df <- data.frame(x = c(cx-s, cx+s, cx+s, cx-s), y = c(cy-s, cy-s, cy+s, cy+s))
    layers[[1]] <- gpoly(sq_df, col)

  } else if (pch == 24L) {
    # Triangle (Deoxyhexose): inscribed in bounding circle R
    # Equilateral triangle with circumradius = R
    a120 <- 2 * pi / 3
    tri_df <- data.frame(
      x = R * cos(c(pi/2, pi/2 - a120, pi/2 + a120)) + cx,
      y = R * sin(c(pi/2, pi/2 - a120, pi/2 + a120)) + cy)
    layers[[1]] <- gpoly(tri_df, col)

  } else if (pch == -1L) {
    # Star (Pentose): outer radius = R, inner = R * 0.52 (plumper)
    star_df <- .dnmb_snfg_star_df(cx, cy, r = R)
    layers[[1]] <- gpoly(star_df, col)

  } else if (pch == -2L) {
    # Crossed square (Hexosamine): diagonal split only, never horizontal
    layers <- .dnmb_snfg_crossed_square_layers(
      cx = cx, cy = cy, s = s, fill_color = col, border_color = bdr, border_lw = blw
    )

  } else if (pch == -3L) {
    # Pentagon (Assigned): circumradius = R
    pent_df <- .dnmb_snfg_pentagon_df(cx, cy, r = R)
    layers[[1]] <- gpoly(pent_df, col)

  } else if (pch == -4L) {
    # Flat rectangle (Di-deoxyhexose): bounding circle = R
    # hw^2 + hh^2 = R^2, aspect ~2:1 → hw = R*cos(atan(0.5)), hh = R*sin(atan(0.5))
    hw <- R * 0.894; hh <- R * 0.447
    rect_df <- data.frame(x = c(cx-hw, cx+hw, cx+hw, cx-hw), y = c(cy-hh, cy-hh, cy+hh, cy+hh))
    layers[[1]] <- gpoly(rect_df, col)

  } else if (pch == -5L) {
    # Flat hexagon (Other): bounding circle = R, squashed vertically
    ang6 <- seq(0, 2 * pi, length.out = 7L)[-(7L)]
    hex_df <- data.frame(x = R * cos(ang6) + cx, y = R * 0.65 * sin(ang6) + cy)
    layers[[1]] <- gpoly(hex_df, col)

  } else if (pch == -6L) {
    # Divided triangle (DeoxyhexNAc): same outline as 24L, vertical split
    a120 <- 2 * pi / 3
    top <- c(cx + R * cos(pi/2), cy + R * sin(pi/2))
    br  <- c(cx + R * cos(pi/2 - a120), cy + R * sin(pi/2 - a120))
    bl  <- c(cx + R * cos(pi/2 + a120), cy + R * sin(pi/2 + a120))
    bot_mid <- c(cx, bl[2])
    left_df  <- data.frame(x = c(top[1], bl[1], bot_mid[1]), y = c(top[2], bl[2], bot_mid[2]))
    right_df <- data.frame(x = c(top[1], bot_mid[1], br[1]), y = c(top[2], bot_mid[2], br[2]))
    outline_df <- data.frame(x = c(top[1], br[1], bl[1]), y = c(top[2], br[2], bl[2]))
    layers[[1]] <- gpoly_na(left_df, "white")
    layers[[2]] <- gpoly_na(right_df, col)
    layers[[3]] <- gpoly(outline_df, fill = NA)
    layers[[4]] <- ggplot2::annotate("segment",
      x = cx, xend = cx, y = top[2], yend = bot_mid[2], color = bdr, linewidth = blw)

  } else if (pch == -7L) {
    # Flat diamond (3,9-dideoxy): bounding circle = R, squashed vertically
    fd_df <- data.frame(
      x = c(cx, cx + R, cx, cx - R),
      y = c(cy + R * 0.5, cy, cy - R * 0.5, cy))
    layers[[1]] <- gpoly(fd_df, col)

  } else if (pch == 23L && is_split) {
    # Divided diamond (Hexuronate): upper triangle colored, lower triangle white
    layers <- .dnmb_snfg_split_diamond_layers(
      cx = cx, cy = cy, r = R, fill_color = col, border_color = bdr, border_lw = blw
    )

  } else if (pch == 23L && !is_split) {
    # Filled diamond (Nonulosonate): same outline as divided diamond
    dm_df <- data.frame(
      x = c(cx, cx + R, cx, cx - R),
      y = c(cy + R, cy, cy - R, cy))
    layers[[1]] <- gpoly(dm_df, col)

  } else {
    layers[[1]] <- ggplot2::annotate("point", x = cx, y = cy,
      shape = pch, size = point_size, fill = col, color = bdr, stroke = blw)
  }
  layers
}

# ---- Render composite SNFG metabolite node ----
# Draws mono/di/trisaccharides with bonds, phospho markers, or plain text
# Returns list of ggplot2 layers
.dnmb_snfg_render_metabolite_node <- function(cx, cy, metabolite, active = TRUE,
                                                symbol_r = 0.33, point_size = 9,
                                                spacing = 0.70) {
  comp <- .dnmb_snfg_metabolite_composition(metabolite)
  layers <- list()
  n_mono <- length(comp$monomers)

  if (n_mono == 0) {
    # Non-sugar metabolite: uniform-size gray circle (same visual weight)
    layers[[1]] <- ggplot2::annotate("point", x = cx, y = cy,
      shape = 21, size = point_size, fill = "#B0BEC5", color = "#CCCCCC", stroke = 0.3)
    return(layers)
  }

  # Center the monomer chain horizontally around cx
  total_w <- (n_mono - 1) * spacing
  x_start <- cx - total_w / 2
  mono_xs <- x_start + (seq_len(n_mono) - 1) * spacing

  # Draw bond lines between monomers — center-to-center (no gap),
  # drawn BEFORE symbols so symbols paint over the line ends
  if (n_mono > 1 && !is.na(comp$linetype)) {
    bond_lty <- comp$linetype
    bond_labels <- comp$bond_labels
    for (j in seq_len(n_mono - 1)) {
      # Determine linetype per bond: for multi-bond (e.g. Raffinose) use per-bond label
      this_lty <- bond_lty
      this_label <- if (length(bond_labels) >= j) bond_labels[j] else bond_labels[1]
      # Alpha = solid, beta = dashed
      if (!is.na(this_label)) {
        this_lty <- if (grepl("^\u03b2", this_label)) "dashed" else "solid"
      }
      bond_df <- data.frame(
        x = mono_xs[j], xend = mono_xs[j + 1],
        y = cy, yend = cy, stringsAsFactors = FALSE)
      layers[[length(layers) + 1]] <- ggplot2::geom_segment(
        data = bond_df,
        ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend),
        color = if (active) "#555555" else "#CCCCCC",
        linewidth = 0.8, linetype = this_lty, inherit.aes = FALSE)
      # Bond label above the line
      if (!is.na(this_label)) {
        mid_x <- (mono_xs[j] + mono_xs[j + 1]) / 2
        layers[[length(layers) + 1]] <- ggplot2::annotate("text",
          x = mid_x, y = cy + symbol_r * 0.9,
          label = this_label, size = 1.5,
          color = if (active) "grey40" else "#CCCCCC",
          fontface = "italic")
      }
    }
  }

  # Draw each monomer SNFG symbol (on top of bond lines)
  for (j in seq_len(n_mono)) {
    sym_layers <- .dnmb_snfg_render_symbol(
      mono_xs[j], cy, comp$monomers[j],
      r = symbol_r, point_size = point_size, active = active)
    for (sl in sym_layers) layers[[length(layers) + 1]] <- sl
  }

  # Phosphorylation marker
  if (comp$phospho) {
    p_x <- mono_xs[n_mono] + symbol_r + 0.08
    layers[[length(layers) + 1]] <- ggplot2::annotate("text",
      x = p_x, y = cy, label = "P",
      size = 2.0, fontface = "bold", color = if (active) "#D84315" else "#CCCCCC")
  }

  layers
}

# ---- Comprehensive metabolism node position table ----
# Glycolysis backbone + PPP branch + all sugar/acid entry points + SNFG annotation
# Layout: compact vertical backbone, sugars fan out left, PPP/pentoses right, acids below
.dnmb_cct_metabolism_nodes <- function() {
  n <- function(id, x, y, label, type = "intermediate",
                sugar_type = "generic", phosphorylated = FALSE) {
    data.frame(id = id, x = x, y = y, label = label, type = type,
               sugar_type = sugar_type, phosphorylated = phosphorylated,
               stringsAsFactors = FALSE)
  }
  # Backbone x-position (center of map)
  bx  <- 7.0
  # PPP x-position (right of backbone)
  px  <- 10.5
  # DHAP offset (left of GA3P)
  dx  <- 5.5

  rbind(
    # ================================================================
    # GLYCOLYSIS BACKBONE — compact vertical, 1.5-unit steps at branch
    # points, 0.8-unit steps where no branches enter
    # ================================================================
    n("Glucose",    bx, 22.0, "Glucose",       "backbone", "glucose"),
    n("Glc-6-P",   bx, 20.5, "Glc-6-P",       "backbone", "glucose",    TRUE),
    n("Fru-6-P",   bx, 19.0, "Fru-6-P",       "backbone", "fructose",   TRUE),
    n("Fru-1,6-BP",bx, 17.8, "Fru-1,6-BP",    "backbone", "fructose",   TRUE),
    n("DHAP",      dx, 16.8, "DHAP",           "backbone", "phospho_sugar", TRUE),
    n("GA3P",      bx, 16.8, "GA3P",           "backbone", "phospho_sugar", TRUE),
    # Lower glycolysis — compact (no external branches)
    n("1,3-BPG",   bx, 15.8, "1,3-BPG",       "backbone", "phospho_sugar", TRUE),
    n("3-PG",      bx, 15.0, "3-PG",           "backbone", "phospho_sugar", TRUE),
    n("2-PG",      bx, 14.2, "2-PG",           "backbone", "phospho_sugar", TRUE),
    n("PEP",       bx, 13.4, "PEP",            "backbone", "phospho_sugar", TRUE),
    # Lower metabolism — spaced for branches
    n("Pyruvate",  bx, 12.0, "Pyruvate",       "backbone", "organic_acid"),
    n("AcCoA",     bx, 10.2, "Acetyl-CoA",     "backbone", "organic_acid"),
    n("TCA",       bx,  8.0, "TCA cycle",      "backbone", "intermediate"),

    # ================================================================
    # PENTOSE PHOSPHATE PATHWAY (PPP) — branch right from Glc-6-P
    # ================================================================
    n("6-PGL",     px, 20.5, "6-PGL",          "ppp", "phospho_sugar", TRUE),
    n("6-PG",      px, 19.8, "6-PG",           "ppp", "phospho_sugar", TRUE),
    n("Ru-5-P",    px, 19.1, "Ru-5-P",         "ppp", "phospho_sugar", TRUE),
    n("R-5-P",     px + 1.2, 18.4, "R-5-P",    "ppp", "ribose",        TRUE),
    n("Xu-5-P",    px, 18.4, "Xu-5-P",         "ppp", "xylose",        TRUE),
    n("S-7-P",     px + 1.2, 17.7, "S-7-P",    "ppp", "phospho_sugar", TRUE),
    n("E-4-P",     px + 1.2, 17.0, "E-4-P",    "ppp", "phospho_sugar", TRUE),

    # ================================================================
    # ENTRY INTERMEDIATES (intracellular conversion steps)
    # Spread vertically to avoid label collision
    # ================================================================
    n("Glc-1-P",    4.5, 21.5, "Glc-1-P",     "entry_intermediate", "glucose",    TRUE),
    n("Gal-1-P",    3.0, 21.0, "Gal-1-P",     "entry_intermediate", "galactose",  TRUE),
    n("GlcNAc-6-P", 4.0, 20.5, "GlcNAc-6-P", "entry_intermediate", "GlcNAc",     TRUE),
    n("GlcN-6-P",   4.0, 19.8, "GlcN-6-P",   "entry_intermediate", "glucosamine", TRUE),
    n("Man-6-P",    4.0, 19.2, "Man-6-P",     "entry_intermediate", "mannose",    TRUE),
    n("Fru-1-P",    4.5, 18.2, "Fru-1-P",     "entry_intermediate", "fructose",   TRUE),
    n("Glycerol-3-P",dx, 15.5, "Glycerol-3-P","entry_intermediate", "glycerol",   TRUE),
    n("KDG",        4.0, 13.8, "KDG",          "entry_intermediate", "organic_acid"),
    n("KDPG",       4.0, 13.0, "KDPG",         "entry_intermediate", "organic_acid"),

    # ================================================================
    # CARBON SOURCES — Glc-6-P entry group (left, upper)
    # ================================================================
    n("glucose_ext",   1.5, 22.0, "Glucose",    "carbon_source", "glucose"),
    n("cellobiose",    0.5, 22.7, "Cellobiose",  "carbon_source", "glucose"),
    n("maltose",       0.5, 21.5, "Maltose",     "carbon_source", "glucose"),
    n("galactose",     0.5, 20.8, "Galactose",   "carbon_source", "galactose"),
    n("trehalose",     0.5, 20.1, "Trehalose",   "carbon_source", "glucose"),
    n("glucose-6-P",   2.5, 20.8, "Glc-6-P*",   "carbon_source", "glucose",     TRUE),

    # ================================================================
    # CARBON SOURCES — Fru-6-P entry group
    # ================================================================
    n("mannose",       1.5, 19.5, "Mannose",     "carbon_source", "mannose"),
    n("NAG",           1.5, 20.2, "NAG",         "carbon_source", "GlcNAc"),
    n("glucosamine",   1.5, 19.8, "Glucosamine", "carbon_source", "glucosamine"),
    n("mannitol",      1.5, 19.0, "Mannitol",    "carbon_source", "mannitol"),
    n("sorbitol",      2.5, 18.6, "Sorbitol",    "carbon_source", "sorbitol"),

    # ================================================================
    # CARBON SOURCES — Fru-1-P / DHAP+GA3P entry
    # ================================================================
    n("fructose",      1.5, 18.2, "Fructose",    "carbon_source", "fructose"),
    n("sucrose",       1.5, 17.5, "Sucrose",     "carbon_source", "glucose"),
    n("lactose",       0.5, 19.4, "Lactose",     "carbon_source", "galactose"),

    # ================================================================
    # CARBON SOURCES — PPP entry (right side)
    # ================================================================
    n("xylose",        px + 3.0, 18.4, "Xylose",     "carbon_source", "xylose"),
    n("xylitol",       px + 3.0, 18.0, "Xylitol",    "carbon_source", "xylitol"),
    n("arabinose",     px + 3.0, 19.1, "Arabinose",  "carbon_source", "arabinose"),
    n("ribose",        px + 3.0, 18.8, "Ribose",     "carbon_source", "ribose"),

    # ================================================================
    # CARBON SOURCES — Nucleoside-derived (far right, near PPP)
    # ================================================================
    n("deoxyinosine",  px + 3.0, 17.6, "dInosine",   "carbon_source", "nucleoside"),
    n("deoxyribose",   px + 3.0, 17.2, "dRibose",    "carbon_source", "nucleoside"),
    n("deoxyribonate", px + 3.0, 16.8, "dRibonate",  "carbon_source", "nucleoside"),
    n("thymidine",     px + 3.0, 16.4, "Thymidine",  "carbon_source", "nucleoside"),

    # ================================================================
    # CARBON SOURCES — GA3P/DHAP entry (deoxy sugars, glycerol, inositol)
    # ================================================================
    n("glycerol",      2.5, 15.5, "Glycerol",     "carbon_source", "glycerol"),
    n("rhamnose",      2.5, 17.2, "Rhamnose",     "carbon_source", "rhamnose"),
    n("fucose",        2.5, 16.8, "Fucose",        "carbon_source", "fucose"),
    n("myoinositol",   2.5, 16.4, "myo-Inositol", "carbon_source", "myoinositol"),

    # ================================================================
    # CARBON SOURCES — Uronic acids / ED pathway (left, mid)
    # ================================================================
    n("gluconate",     1.5, 14.2, "Gluconate",     "carbon_source", "GlcA"),
    n("glucuronate",   1.5, 13.6, "Glucuronate",   "carbon_source", "GlcA"),
    n("galacturonate", 1.5, 13.0, "Galacturonate", "carbon_source", "GalA"),

    # ================================================================
    # CARBON SOURCES — Pyruvate entry (left, below PEP)
    # ================================================================
    n("L-lactate",     4.5, 12.5, "L-Lactate",    "carbon_source", "organic_acid"),
    n("D-lactate",     4.5, 12.0, "D-Lactate",    "carbon_source", "organic_acid"),
    n("pyruvate_ext",  4.5, 11.5, "Pyruvate*",    "carbon_source", "organic_acid"),

    # ================================================================
    # CARBON SOURCES — AcCoA entry
    # ================================================================
    n("acetate",       4.5, 10.6, "Acetate",       "carbon_source", "organic_acid"),
    n("ethanol",       4.5, 10.2, "Ethanol",       "carbon_source", "organic_acid"),
    n("propionate",    4.5, 9.8,  "Propionate",    "carbon_source", "organic_acid"),

    # ================================================================
    # CARBON SOURCES — TCA entry (left, bottom)
    # ================================================================
    n("citrate",       4.5, 8.8, "Citrate",        "carbon_source", "organic_acid"),
    n("succinate",     4.5, 8.4, "Succinate",       "carbon_source", "organic_acid"),
    n("fumarate",      4.5, 8.0, "Fumarate",       "carbon_source", "organic_acid"),
    n("L-malate",      4.5, 7.6, "L-Malate",       "carbon_source", "organic_acid"),
    n("2-oxoglutarate",4.5, 7.2, "2-OG",           "carbon_source", "organic_acid"),

    # ================================================================
    # CARBON SOURCES — Amino acids -> Pyruvate (right of backbone)
    # ================================================================
    n("alanine",       9.0, 12.4, "Alanine",       "carbon_source", "amino_acid"),
    n("serine",        9.0, 12.0, "Serine",        "carbon_source", "amino_acid"),
    n("D-alanine",     9.0, 11.6, "D-Alanine",     "carbon_source", "amino_acid"),
    n("D-serine",      9.0, 11.2, "D-Serine",      "carbon_source", "amino_acid"),
    n("threonine",     9.0, 10.8, "Threonine",     "carbon_source", "amino_acid"),

    # ================================================================
    # CARBON SOURCES — Amino acids -> AcCoA (right, mid)
    # ================================================================
    n("lysine",        10.5, 10.8, "Lysine",       "carbon_source", "amino_acid"),
    n("leucine",       10.5, 10.4, "Leucine",      "carbon_source", "amino_acid"),
    n("isoleucine",    10.5, 10.0, "Isoleucine",   "carbon_source", "amino_acid"),
    n("valine",        10.5, 9.6, "Valine",        "carbon_source", "amino_acid"),
    n("phenylalanine", 11.8, 10.8, "Phe",          "carbon_source", "amino_acid"),
    n("tyrosine",      11.8, 10.4, "Tyr",          "carbon_source", "amino_acid"),
    n("tryptophan",    11.8, 10.0, "Trp",          "carbon_source", "amino_acid"),

    # ================================================================
    # CARBON SOURCES — Amino acids -> TCA (right, bottom)
    # ================================================================
    n("glutamate",     9.0, 8.8,  "Glutamate",     "carbon_source", "amino_acid"),
    n("aspartate",     9.0, 8.4,  "Aspartate",     "carbon_source", "amino_acid"),
    n("asparagine",    9.0, 8.0,  "Asparagine",    "carbon_source", "amino_acid"),
    n("proline",       9.0, 7.6,  "Proline",       "carbon_source", "amino_acid"),
    n("arginine",      10.5, 8.4, "Arginine",      "carbon_source", "amino_acid"),
    n("histidine",     10.5, 8.0, "Histidine",     "carbon_source", "amino_acid"),
    n("citrulline",    10.5, 7.6, "Citrulline",    "carbon_source", "amino_acid"),

    # ================================================================
    # CARBON SOURCES — Aromatics (far right, mid)
    # ================================================================
    n("4-hydroxybenzoate", 11.8, 9.6, "4-HBA",        "carbon_source", "aromatic"),
    n("phenylacetate",     11.8, 9.2, "Phenylacetate", "carbon_source", "aromatic"),

    # ================================================================
    # CAZy bridge nodes — compact row at top-left, linked to disaccharides
    # ================================================================
    n("CAZy_GH", -0.8, 23.0, "GH",  "cazy_bridge", "generic"),
    n("CAZy_CE",  0.5, 23.0, "CE",  "cazy_bridge", "generic"),
    n("CAZy_GT",  1.8, 23.0, "GT",  "cazy_bridge", "generic"),
    n("CAZy_AA",  3.1, 23.0, "AA",  "cazy_bridge", "generic")
  )
}

# SNFG-annotated node data: merges position table with SNFG shape/color/fill
.dnmb_cct_nodes_styled <- function() {
  nodes <- .dnmb_cct_metabolism_nodes()
  snfg  <- .dnmb_snfg_lookup(nodes$sugar_type)
  nodes$snfg_shape <- snfg$snfg_shape
  nodes$snfg_color <- snfg$snfg_color
  nodes$snfg_fill  <- snfg$snfg_fill
  # Phosphorylated sugars get a 'P' label suffix indicator
  nodes$p_label <- ifelse(nodes$phosphorylated, "P", "")
  nodes
}

# Metabolism backbone edges (glycolysis + PPP internal connections)
.dnmb_cct_backbone_edges <- function(nodes) {
  nx <- stats::setNames(nodes$x, nodes$id)
  ny <- stats::setNames(nodes$y, nodes$id)
  e <- function(from, to, pathway = "glycolysis") {
    if (!from %in% names(nx) || !to %in% names(nx)) return(NULL)
    data.frame(from = from, to = to,
               from_x = nx[from], from_y = ny[from],
               to_x = nx[to], to_y = ny[to],
               pathway = pathway, stringsAsFactors = FALSE)
  }
  do.call(rbind, Filter(Negate(is.null), list(
    # --- Glycolysis ---
    e("Glucose",     "Glc-6-P"),
    e("Glc-6-P",    "Fru-6-P"),
    e("Fru-6-P",    "Fru-1,6-BP"),
    e("Fru-1,6-BP", "GA3P"),
    e("Fru-1,6-BP", "DHAP"),
    e("DHAP",       "GA3P"),
    e("GA3P",       "1,3-BPG"),
    e("1,3-BPG",    "3-PG"),
    e("3-PG",       "2-PG"),
    e("2-PG",       "PEP"),
    e("PEP",        "Pyruvate"),
    e("Pyruvate",   "AcCoA"),
    e("AcCoA",      "TCA"),
    # --- PPP oxidative ---
    e("Glc-6-P",    "6-PGL",    "ppp"),
    e("6-PGL",      "6-PG",     "ppp"),
    e("6-PG",       "Ru-5-P",   "ppp"),
    e("Ru-5-P",     "R-5-P",    "ppp"),
    e("Ru-5-P",     "Xu-5-P",   "ppp"),
    # --- PPP non-oxidative ---
    e("R-5-P",      "S-7-P",    "ppp"),
    e("Xu-5-P",     "S-7-P",    "ppp"),
    e("S-7-P",      "E-4-P",    "ppp"),
    e("Xu-5-P",     "GA3P",     "ppp"),
    e("S-7-P",      "Fru-6-P",  "ppp"),
    e("E-4-P",      "Fru-6-P",  "ppp"),
    # --- Entner-Doudoroff ---
    e("KDG",        "KDPG",     "ed"),
    e("KDPG",       "Pyruvate", "ed"),
    e("KDPG",       "GA3P",     "ed")
  )))
}

# Entry connections: carbon source -> first intracellular intermediate -> backbone
.dnmb_cct_entry_edges <- function(nodes) {
  nx <- stats::setNames(nodes$x, nodes$id)
  ny <- stats::setNames(nodes$y, nodes$id)
  e <- function(from, to) {
    if (!from %in% names(nx) || !to %in% names(nx)) return(NULL)
    data.frame(from = from, to = to,
               from_x = nx[from], from_y = ny[from],
               to_x = nx[to], to_y = ny[to],
               stringsAsFactors = FALSE)
  }
  do.call(rbind, Filter(Negate(is.null), list(
    # --- Glc-6-P entry ---
    e("glucose_ext", "Glucose"),
    e("maltose",     "Glc-1-P"),   e("Glc-1-P", "Glc-6-P"),
    e("galactose",   "Gal-1-P"),   e("Gal-1-P", "Glc-1-P"),
    e("trehalose",   "Glucose"),
    e("cellobiose",  "Glucose"),
    e("glucose-6-P", "Glc-6-P"),
    # --- Fru-6-P entry ---
    e("mannose",     "Man-6-P"),   e("Man-6-P",     "Fru-6-P"),
    e("NAG",         "GlcNAc-6-P"),e("GlcNAc-6-P", "GlcN-6-P"),
    e("GlcN-6-P",   "Fru-6-P"),
    e("glucosamine", "GlcN-6-P"),
    e("mannitol",    "Fru-6-P"),
    # --- Fru-1-P entry ---
    e("fructose",    "Fru-1-P"),   e("Fru-1-P", "DHAP"),
    e("Fru-1-P",    "GA3P"),
    e("sorbitol",    "Fru-6-P"),
    # --- Sucrose splits ---
    e("sucrose",     "Glucose"),
    e("sucrose",     "Fru-1-P"),
    # --- Lactose ---
    e("lactose",     "Glucose"),
    e("lactose",     "Gal-1-P"),
    # --- PPP entry ---
    e("xylose",      "Xu-5-P"),
    e("xylitol",     "Xu-5-P"),
    e("arabinose",   "Ru-5-P"),
    e("ribose",      "R-5-P"),
    # --- Nucleoside-derived -> PPP / GA3P ---
    e("deoxyinosine",  "R-5-P"),
    e("deoxyribose",   "R-5-P"),
    e("deoxyribonate", "GA3P"),
    e("thymidine",     "R-5-P"),
    # --- Sugar alcohols ---
    e("glycerol",      "Glycerol-3-P"), e("Glycerol-3-P", "DHAP"),
    e("myoinositol",   "GA3P"),
    # --- Deoxy sugars ---
    e("rhamnose",    "DHAP"),
    e("fucose",      "DHAP"),
    # --- Uronic acids (Entner-Doudoroff variant) ---
    e("gluconate",     "KDG"),
    e("glucuronate",   "KDG"),
    e("galacturonate", "KDG"),
    # --- Pyruvate/TCA entry ---
    e("L-lactate",     "Pyruvate"),
    e("D-lactate",     "Pyruvate"),
    e("pyruvate_ext",  "Pyruvate"),
    e("acetate",       "AcCoA"),
    e("ethanol",       "AcCoA"),
    e("propionate",    "AcCoA"),
    e("citrate",       "TCA"),
    e("succinate",     "TCA"),
    e("fumarate",      "TCA"),
    e("L-malate",      "TCA"),
    e("2-oxoglutarate","TCA"),
    # --- Amino acids ---
    e("glutamate",     "TCA"),
    e("aspartate",     "TCA"),
    e("asparagine",    "TCA"),
    e("proline",       "TCA"),
    e("arginine",      "TCA"),
    e("histidine",     "TCA"),
    e("citrulline",    "TCA"),
    e("alanine",       "Pyruvate"),
    e("serine",        "Pyruvate"),
    e("D-alanine",     "Pyruvate"),
    e("D-serine",      "Pyruvate"),
    e("threonine",     "Pyruvate"),
    e("lysine",        "AcCoA"),
    e("leucine",       "AcCoA"),
    e("isoleucine",    "AcCoA"),
    e("valine",        "AcCoA"),
    e("phenylalanine", "AcCoA"),
    e("tyrosine",      "AcCoA"),
    e("tryptophan",    "AcCoA"),
    # --- Aromatics ---
    e("4-hydroxybenzoate", "AcCoA"),
    e("phenylacetate",     "AcCoA")
  )))
}

# Legacy wrapper — existing code calls .dnmb_gapmind_carbon_cazy_nodes()
# Maps new comprehensive node table back to the old schema for compatibility
#' Internal SNFG symbol drawing helpers
#'
#' Low-level geometry and drawing helpers used for SNFG-style carbohydrate
#' symbols in DNMB figures.
#'
#' @name dnmb_internal_snfg_symbols
#' @keywords internal
#' @noRd
NULL
