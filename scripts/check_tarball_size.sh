#!/usr/bin/env sh
# Fail if the R source tarball exceeds CRAN's incoming size ceiling.
#
# Why this exists: the budget is not a nice-to-have. CRAN's incoming check flags any
# tarball over 5,000,000 bytes (`_R_CHECK_CRAN_INCOMING_TARBALL_THRESHOLD_`, which
# defaults to `5e6`), and 5 MB is the documented policy ceiling for a source package.
# This package lives at ~4.95 MB, i.e. about 45 KB -- under one percent -- below that
# line, because it ships both a 2.9 MB database export tree and a 3.6 MB pinned upstream
# tree that every differential gate diffs against. Any future addition to `R/`, `man/`,
# `inst/db` or `vignettes/` can push it over again, and nothing else in the repository
# fails when that happens: `R CMD check` reports the size as a NOTE, which is easy to
# miss in a 700-line log. So this script fails loudly, and CI runs it right after
# `R CMD build`, before the hours-long check.
#
# What to do when it fires: do not just raise the ceiling. The margin is the point.
# Shrink the contents instead -- see `.Rbuildignore`, which documents what is already
# excluded and why. The next candidates, measured, are `inst/db/*/symbols.txt` (~374 KB,
# md5-verified by `r-core/src/db.rs`, so the loader has to agree) and the vendored crate
# archive (`src/rust/r-build/vendor.tar.xz`, 962 KB; cargo checksums every vendored file,
# so pruning is a `vendor_rust_deps.sh` change, not a `tar` one).
set -eu

LIMIT=5000000
TARBALL="${1:-}"
if [ -z "$TARBALL" ]; then
  # shellcheck disable=SC2012
  TARBALL="$(ls -t CellChat_*.tar.gz 2>/dev/null | head -n 1 || true)"
fi
if [ -z "$TARBALL" ] || [ ! -f "$TARBALL" ]; then
  echo "check_tarball_size: no source tarball found (expected CellChat_*.tar.gz)" >&2
  exit 1
fi

size="$(wc -c < "$TARBALL" | tr -d '[:space:]')"
echo "check_tarball_size: $TARBALL is $size bytes against a $LIMIT-byte ceiling"
if [ "$size" -gt "$LIMIT" ]; then
  echo "check_tarball_size: OVER the CRAN incoming size ceiling by $((size - LIMIT)) bytes." >&2
  echo "check_tarball_size: largest contents:" >&2
  tar -tzvf "$TARBALL" 2>/dev/null \
    | awk '{print $3, $6}' | sort -rn | awk 'NR<=10 {print "  " $0}' >&2
  exit 1
fi
echo "check_tarball_size: under by $((LIMIT - size)) bytes"
