suppressPackageStartupMessages(library(CellChat))

stopifnot(cellchatrs_available())
stopifnot(methods::isClass("CellChat"))

## The installed package must deserialize CellChat objects using its registered S4 class and
## provide the pinned fallback without changing the caller's search path or global workspace.
local({
  global_before <- ls(.GlobalEnv, all.names = TRUE)
  search_before <- search()
  warnings <- character()
  cached <- get("cellchatrs_upstream_cached", envir = asNamespace("CellChat"))
  reference <- withCallingHandlers(
    cached(),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  stopifnot(
    is.environment(reference),
    is.function(get("computeCommunProb", envir = reference)),
    identical(search_before, search()),
    identical(global_before, ls(.GlobalEnv, all.names = TRUE)),
    length(warnings) == 0L
  )

  object <- methods::new("CellChat")
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path)
  stopifnot(identical(object, readRDS(path)))
})
