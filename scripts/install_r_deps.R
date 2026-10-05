#!/usr/bin/env Rscript
# Install R packages into an explicit directory, without reading any environment
# variable for the location.
#
# Why this exists instead of `lib = Sys.getenv("R_LIBS_USER")`: something in the CI
# R startup (an `Renviron.site` or `Rprofile.site` from the toolchain setup, most
# likely) appends `:/opt/R/<ver>/lib/R/library` to `R_LIBS_USER` before any package
# code runs, so `Sys.getenv()` inside `Rscript -e` returns a colon-separated pair
# even when the shell exported a single path. Every installer then reports the
# combined string as not writable and installs nothing, while fixture generation --
# which reads the same variable -- looks in directories that were never written.
# Passing the directory as a command-line argument sidesteps all of it: no
# environment variable is consulted for the location at any point.
#
# Usage:
#   Rscript scripts/install_r_deps.R --lib /abs/path --cran pkg1,pkg2 --bioc pkg3
#   Rscript scripts/install_r_deps.R --lib /abs/path --pin collapse=2.1.8,dplyr=1.2.1
#
# `--cran` installs latest via `install.packages` with recursive Depends/Imports, so
# transitive hard dependencies arrive too (base R resolves the closure; `remotes`
# reported success while leaving `sass`, `rstatix` and `htmlwidgets` uninstalled).
# `--bioc` goes through BiocManager. `--pin` uses `remotes::install_version` for the
# bit-identity-critical packages whose exact output the goldens pin; it needs
# `remotes`, installed on first use. All three flags may be combined; `--lib` is
# created if missing.
args <- commandArgs(trailingOnly = TRUE)
get_flag <- function(name) {
  # Accept `--name=value` and `--name value` alike; CI calls use a space because the
  # value is a shell variable, and a missing value is a hard error rather than a
  # silent install-everything.
  eq <- grep(paste0("^", name, "="), args, value = TRUE)
  if (length(eq)) return(sub(paste0("^", name, "="), "", eq[1]))
  i <- match(name, args)
  if (!is.na(i) && i < length(args)) return(args[i + 1])
  NULL
}
lib <- get_flag("--lib")
if (is.null(lib) || !nzchar(lib)) stop("--lib PATH is required", call. = FALSE)
lib <- normalizePath(lib, mustWork = FALSE)
dir.create(lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

cran <- get_flag("--cran")
if (!is.null(cran) && nzchar(cran)) {
  pkgs <- strsplit(cran, ",")[[1]]
  utils::install.packages(pkgs, repos = "https://cloud.r-project.org",
                          lib = lib,
                          dependencies = c("Depends", "Imports", "LinkingTo"))
}

bioc <- get_flag("--bioc")
if (!is.null(bioc) && nzchar(bioc)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    utils::install.packages("BiocManager", repos = "https://cloud.r-project.org", lib = lib)
  BiocManager::install(strsplit(bioc, ",")[[1]], ask = FALSE, update = FALSE, lib = lib)
}

pin <- get_flag("--pin")
if (!is.null(pin) && nzchar(pin)) {
  if (!requireNamespace("remotes", quietly = TRUE))
    utils::install.packages("remotes", repos = "https://cloud.r-project.org", lib = lib)
  for (spec in strsplit(pin, ",")[[1]]) {
    kv <- strsplit(spec, "=")[[1]]
    remotes::install_version(kv[1], version = kv[2], lib = lib, upgrade = "never")
  }
}

cat("install_r_deps: library", lib, "now holds", length(rownames(utils::installed.packages(lib.loc = lib))), "packages\n")
