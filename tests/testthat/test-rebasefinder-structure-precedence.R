test_that("user and ESMFold structures outrank managed ProMod3 models", {
  output <- tempfile("rebasefinder-output-")
  validation_root <- tempfile("rebasefinder-validation-")
  user_dir <- file.path(output, "query_structures")
  managed_dir <- file.path(
    output, "dnmb_module_rebasefinder", "promod3_query_structures"
  )
  validation_user_dir <- file.path(validation_root, "query_structures")
  esmfold_dir <- file.path(validation_root, "esmfold_query_structures")
  dir.create(user_dir, recursive = TRUE)
  dir.create(managed_dir, recursive = TRUE)
  dir.create(validation_user_dir, recursive = TRUE)
  dir.create(esmfold_dir, recursive = TRUE)
  on.exit(unlink(c(output, validation_root), recursive = TRUE, force = TRUE), add = TRUE)

  validation_path <- file.path(validation_root, "foldseek.tsv")
  writeLines("query\thit", validation_path)
  dirs <- normalizePath(DNMB:::.dnmb_rebasefinder_structure_dirs(
    output, validation_path
  ))

  independent <- normalizePath(c(user_dir, validation_user_dir, esmfold_dir))
  managed <- normalizePath(managed_dir)
  expect_true(all(match(independent, dirs) < match(managed, dirs)))

  user_structure <- file.path(user_dir, "user_query.cif")
  validation_user_structure <- file.path(validation_user_dir, "validation_query.cif")
  esmfold_structure <- file.path(esmfold_dir, "esm_query.cif")
  managed_structures <- file.path(
    managed_dir,
    c("user_query.pdb", "validation_query.pdb", "esm_query.pdb")
  )
  writeLines("independent", user_structure)
  writeLines("independent", validation_user_structure)
  writeLines("independent", esmfold_structure)
  invisible(base::lapply(managed_structures, function(path) writeLines("managed", path)))

  expect_identical(
    DNMB:::.dnmb_rebasefinder_structure_file_for_query("user_query", dirs),
    normalizePath(user_structure)
  )
  expect_identical(
    DNMB:::.dnmb_rebasefinder_structure_file_for_query("validation_query", dirs),
    normalizePath(validation_user_structure)
  )
  expect_identical(
    DNMB:::.dnmb_rebasefinder_structure_file_for_query("esm_query", dirs),
    normalizePath(esmfold_structure)
  )
})
