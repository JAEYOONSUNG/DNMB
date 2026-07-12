test_that("cytoplasmic sugar nodes use compact labels inside their glyphs", {
  inside <- DNMB:::.dnmb_cct_cytoplasmic_inside_label

  expect_identical(inside("Glucose", "Glucose", "backbone", "glucose"), "Glc")
  expect_identical(
    inside("GlcNAc-6-P", "GlcNAc-6-P", "entry_intermediate", "glcnac"),
    "GlcNAc-\n6-P"
  )
  expect_identical(
    inside("Fru-1,6-BP", "Fru-1,6-BP", "backbone", "fructose"),
    "Fru-\n1,6-BP"
  )
  expect_identical(inside("R-5-P", "R-5-P", "ppp", "ribose"), "R-5-P")
  expect_identical(inside("Glycerol-3-P", "Glycerol-3-P", "entry_intermediate", "glycerol"), "Gro-3-P")
  expect_identical(inside("Pyruvate", "Pyruvate", "backbone", "intermediate"), "")
  expect_identical(inside("Citrate", "Citrate", "tca", "intermediate"), "")
})

test_that("SNFG abbreviations use unambiguous standard residue names", {
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("NAG"), "GlcNAc")
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("glucosamine"), "GlcN")
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("glca"), "GlcA")
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("gala"), "GalA")
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("deoxyribose"), "dRib")
  expect_identical(DNMB:::.dnmb_snfg_abbreviation_v2("glycerol"), "Gro")
})

test_that("carbon-source glyphs contain abbreviations instead of full names", {
  layer_labels <- function(plot) {
    unlist(lapply(plot$layers, function(layer) {
      if (!is.data.frame(layer$data) || !"label" %in% names(layer$data)) return(character(0))
      as.character(layer$data$label)
    }), use.names = FALSE)
  }

  maltose <- DNMB:::.dnmb_cct_render_carbon_source_node(
    ggplot2::ggplot(), 1, 1, "Maltose"
  )
  nag <- DNMB:::.dnmb_cct_render_carbon_source_node(
    ggplot2::ggplot(), 1, 1, "NAG"
  )

  expect_true("Mal" %in% layer_labels(maltose))
  expect_false("Maltose" %in% layer_labels(maltose))
  expect_true("Glc\nNAc" %in% layer_labels(nag))
  expect_false("NAG" %in% layer_labels(nag))
})

test_that("SNFG legend is generated from every rendered sugar type", {
  nodes <- data.frame(
    sugar_type = c("glucose", "intermediate", "deoxyribose", "glucosamine", "glycerol"),
    stringsAsFactors = FALSE
  )
  branches <- list(
    Xylan = list(list(sugar = "arabinose"))
  )
  used <- DNMB:::.dnmb_cct_used_snfg_types(
    cyto_nodes = nodes,
    source_ids = c("Mannitol", "Lactose"),
    monomer_ids = "gala",
    extra_chains = list(Xylan = c("xylose", "xylose")),
    extra_branches = branches
  )
  legend <- DNMB:::.dnmb_cct_snfg_legend_data(used)

  expected <- c(
    "glucose", "deoxyribose", "glucosamine", "glycerol", "mannitol",
    "galactose", "gala", "xylose", "arabinose"
  )
  expect_setequal(used, expected)
  expect_setequal(legend$sugar, used)
  expect_false("intermediate" %in% legend$sugar)
  expect_identical(
    legend$label[legend$sugar == "deoxyribose"],
    "2-Deoxy-D-ribose (dRib)"
  )
  expect_identical(
    legend$label[legend$sugar == "glucosamine"],
    "D-Glucosamine (GlcN)"
  )
  expect_identical(
    legend$label[legend$sugar == "mannitol"],
    "D-Mannitol (Mtl)"
  )
})
