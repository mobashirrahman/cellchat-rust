#!/usr/bin/env bash
# Exercise the distribution boundary outside the checkout, with a fresh build target.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
check_dir="$(mktemp -d "${TMPDIR:-/tmp}/cellchat-source-check.XXXXXX")"
trap 'rm -rf "$check_dir"' EXIT
mkdir -p "$check_dir/library"
(
  cd "$check_dir"
  R CMD build --no-build-vignettes "$root"
  CARGO_TARGET_DIR="$check_dir/target" CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-2}" \
    R CMD INSTALL --library="$check_dir/library" ./*.tar.gz
)
R_LIBS="$check_dir/library${R_LIBS:+:$R_LIBS}" \
  Rscript "$root/tests/test-package.R"
