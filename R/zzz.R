## @useDynLib cellchatrs, .registration = TRUE

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  dll <- getLoadedDLLs()[["cellchatrs"]]

  ## extendr 0.9 exposes compiled entry points as registered `.Call` routines, not as R
  ## functions. Generate thin R wrappers from the module metadata and install them into
  ## the package namespace, which is what `extendr-engine` does at build time.
  ## Doing it here (rather than a generated R/ file) means the wrappers can never drift
  ## out of sync with the Rust module list.
  gen <- getNativeSymbolInfo("wrap__make_cellchatrs_wrappers", PACKAGE = dll)
  ## use_symbols = FALSE -> wrappers call .Call("sym", PACKAGE = "cellchatrs"),
  ## which resolves through useDynLib(.registration = TRUE) without needing
  ## NativeSymbolInfo objects in the namespace.
  src <- .Call(gen$address, FALSE, "cellchatrs")
  for (nm in ls(ns, all.names = TRUE)) {
    if (startsWith(nm, "wrap__") || nm %in% c("R_init_r_bindings", "R_init_libr_bindings")) {
      try(rm(list = nm, envir = ns), silent = TRUE)
    }
  }
  eval(parse(text = src), envir = ns)

  ## rayon cannot read R's environment, so size the pool from here. This only has an
  ## effect before the first kernel call: rayon global pools are built once per process.
  ## NOTE: parallel::detectCores() reports 16 for logical = FALSE on this Zen3 host
  ## (8c/16t), so it must not be trusted for benchmarking. Always pin explicitly with
  ## CELLCHATRS_THREADS, or run under taskset.
  n <- suppressWarnings(as.integer(Sys.getenv("CELLCHATRS_THREADS", "0")))
  if (is.na(n) || n <= 0L) n <- parallel::detectCores(logical = TRUE) %||% 1L
  ## `request_num_threads` only records the request; building the pool here would
  ## make every later set_num_threads() a no-op.
  try(get("request_num_threads", envir = ns)(n), silent = TRUE)
  invisible()
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Is the Rust kernel loaded and callable?
#' @export
cellchatrs_available <- function() {
  is.function(get0("computeCommunProb", envir = asNamespace("cellchatrs"), inherits = FALSE))
}

#' Number of worker threads the Rust kernel is using.
#' @export
cellchatrs_threads <- function() {
  get("get_num_threads", envir = asNamespace("cellchatrs"))()
}

