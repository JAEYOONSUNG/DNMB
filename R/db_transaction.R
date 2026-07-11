.dnmb_transactional_replace <- function(staged_paths,
                                        destination_paths,
                                        retire_paths = character(),
                                        rename_file = base::file.rename) {
  staged_paths <- base::as.character(staged_paths)
  destination_paths <- base::as.character(destination_paths)
  retire_paths <- base::setdiff(base::as.character(retire_paths), destination_paths)

  if (!base::length(staged_paths) ||
      base::length(staged_paths) != base::length(destination_paths) ||
      base::anyDuplicated(destination_paths) > 0L) {
    return(list(ok = FALSE, detail = "Invalid staged generation mapping."))
  }
  missing_staged <- staged_paths[!base::file.exists(staged_paths)]
  if (base::length(missing_staged)) {
    return(list(
      ok = FALSE,
      detail = base::paste0("Staged generation is incomplete: ", missing_staged[[1]])
    ))
  }

  existing_paths <- base::unique(base::c(
    destination_paths[base::file.exists(destination_paths)],
    retire_paths[base::file.exists(retire_paths)]
  ))
  backup_paths <- stats::setNames(
    base::vapply(existing_paths, function(path) {
      base::tempfile(
        pattern = base::paste0(".", base::basename(path), ".dnmb-backup-"),
        tmpdir = base::dirname(path)
      )
    }, character(1)),
    existing_paths
  )
  moved_to_backup <- character()
  installed <- character()

  restore_generation <- function() {
    if (base::length(installed)) {
      base::unlink(installed, recursive = TRUE, force = TRUE)
    }
    restore_ok <- TRUE
    for (path in base::rev(moved_to_backup)) {
      backup <- backup_paths[[path]]
      if (base::file.exists(backup)) {
        restored <- base::suppressWarnings(base::file.rename(backup, path))
        restore_ok <- base::isTRUE(restored) && restore_ok
      }
    }
    restore_ok
  }

  for (path in existing_paths) {
    moved <- tryCatch(
      base::suppressWarnings(rename_file(path, backup_paths[[path]])),
      error = function(e) FALSE
    )
    if (!base::isTRUE(moved)) {
      restored <- restore_generation()
      return(list(
        ok = FALSE,
        detail = base::paste0(
          "Could not stage the previous generation for replacement: ", path,
          if (!restored) "; rollback was incomplete" else ""
        )
      ))
    }
    moved_to_backup <- base::c(moved_to_backup, path)
  }

  for (i in base::seq_along(destination_paths)) {
    destination <- destination_paths[[i]]
    base::dir.create(base::dirname(destination), recursive = TRUE, showWarnings = FALSE)
    moved <- tryCatch(
      base::suppressWarnings(rename_file(staged_paths[[i]], destination)),
      error = function(e) FALSE
    )
    if (!base::isTRUE(moved)) {
      restored <- restore_generation()
      return(list(
        ok = FALSE,
        detail = base::paste0(
          "Could not commit the staged generation to ", destination,
          if (!restored) "; rollback was incomplete" else "; the previous generation was restored"
        )
      ))
    }
    installed <- base::c(installed, destination)
  }

  if (base::length(backup_paths)) {
    base::unlink(base::unname(backup_paths), recursive = TRUE, force = TRUE)
  }
  list(ok = TRUE, detail = "Staged generation committed.")
}
