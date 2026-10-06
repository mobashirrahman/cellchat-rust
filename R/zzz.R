## @useDynLib CellChat, .registration = TRUE

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  dll <- getLoadedDLLs()[["CellChat"]]

  ## extendr 0.9 exposes compiled entry points as registered `.Call` routines, not as R
  ## functions. Generate thin R wrappers from the module metadata and install them into
  ## the package namespace, which is what `extendr-engine` does at build time.
  ## Doing it here (rather than a generated R/ file) means the wrappers can never drift
  ## out of sync with the Rust module list.
  gen <- getNativeSymbolInfo("wrap__make_cellchatrs_wrappers", PACKAGE = dll)
  ## use_symbols = FALSE -> wrappers call .Call("sym", PACKAGE = "CellChat"),
  ## which resolves through useDynLib(.registration = TRUE) without needing
  ## NativeSymbolInfo objects in the namespace.
  src <- .Call(gen$address, FALSE, "CellChat")
  for (nm in ls(ns, all.names = TRUE)) {
    if (startsWith(nm, "wrap__") || nm %in% c("R_init_r_bindings", "R_init_libr_bindings")) {
      try(rm(list = nm, envir = ns), silent = TRUE)
    }
  }
  eval(parse(text = src), envir = ns)

  ## rayon cannot read R's environment, so size the pool from here. This only has an
  ## effect before the first kernel call: rayon global pools are built once per process.
  ## Containers and shared runners often report the host's full CPU count rather than the
  ## allocation available to this process. Keep the default modest; users can opt into a larger
  ## pool explicitly with CELLCHATRS_THREADS.
  n <- suppressWarnings(as.integer(Sys.getenv("CELLCHATRS_THREADS", "0")))
  if (is.na(n) || n <= 0L) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores) || cores < 1L) cores <- 1L
    n <- min(2L, cores)
  }
  ## `request_num_threads` only records the request; building the pool here would
  ## make every later set_num_threads() a no-op.
  try(get("request_num_threads", envir = ns)(n), silent = TRUE)

  ## Upstream declares `LazyData: TRUE`, so its tutorials write `CellChatDB.human` bare. This
  ## package cannot use LazyData: the five datasets ship once, inside `inst/upstream`, because a
  ## second copy in `data/` is what put the source tarball over CRAN's 5 MB limit, and
  ## `.Rbuildignore` drops `data/` from every built install. Bind the names here instead, as
  ## promises that read the bundled `.rda` on first use -- the same thing LazyData does, from the
  ## copy that is already shipped. Nothing is read at load time. `data(CellChatDB.human)` still
  ## warns in a built install, since `data()` only looks in a `data/` directory; the object is
  ## there regardless.
  bundled <- file.path(libname, pkgname, "upstream", "CellChat-75253cd0", "data")
  for (nm in .cellchatrs_datasets) {
    local({
      name <- nm
      delayedAssign(name, {
        e <- new.env(parent = emptyenv())
        load(file.path(bundled, paste0(name, ".rda")), envir = e)
        get(name, envir = e, inherits = FALSE)
      }, assign.env = ns)
    })
  }
  invisible()
}

.cellchatrs_datasets <- c("CellChatDB.human", "CellChatDB.mouse", "CellChatDB.zebrafish",
                          "PPI.human", "PPI.mouse")

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Is the Rust kernel loaded and callable?
#' @export
cellchatrs_available <- function() {
  is.function(get0("computeCommunProb", envir = asNamespace("CellChat"), inherits = FALSE))
}

#' Number of worker threads the Rust kernel is using.
#' @export
cellchatrs_threads <- function() {
  get("get_num_threads", envir = asNamespace("CellChat"))()
}
