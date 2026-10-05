#!/bin/sh
# Regenerate the golden corpora from the pinned upstream CellChat tree.
#
# Three of these are large (rng_stress.txt ~24 MB, f80_ref.txt ~19 MB, stats_vectors.bin
# ~9.5 MB) and are excluded from the source tarball by .Rbuildignore to stay under CRAN's
# 5 MB limit, so a fresh clone needs this script before `cargo test` will pass.
#
# Usage:  sh scripts/gen_fixtures.sh
set -e
cd "$(dirname "$0")/.."
CC="${CELLCHAT_SRC:-../CellChat}"
if [ ! -d "$CC/R" ]; then
  echo "pinned upstream not found at $CC; set CELLCHAT_SRC" >&2
  exit 1
fi
# `R_LIBS_USER`, not `R_LIBS`: the latter is not an R variable at all, so the old
# `R_LIBS=.rlib` prefix was silently ignored and every generator ran with the default
# library paths. Search the ambient value, CI's dedicated package dir, and the local
# `.rlib`, in that order -- `r-lib/actions/setup-r` overrides `R_LIBS_USER` to its own
# temp library while installs may have targeted `$R_PKG_LIB`, so no single one of these
# is guaranteed to hold the packages. R searches a colon-separated `R_LIBS_USER` left
# to right, then the system libraries. Set once: re-appending inside the loop would
# grow the path on every generator.
R_LIBS_USER="${R_LIBS_USER:-}:${R_PKG_LIB:-}:$(pwd)/.rlib"
# Drop empty entries so the intent is explicit (a leading or doubled colon is an empty
# first entry, which R would otherwise resolve unpredictably).
R_LIBS_USER="$(printf '%s' "$R_LIBS_USER" | tr ':' '\n' | grep -v '^$' | paste -sd: -)"
export R_LIBS_USER
echo "gen_fixtures: R_LIBS_USER=$R_LIBS_USER"
for g in tests/parity/gen_*.R; do
  # `gen_tutorial_netp_slice.R` is documented one-off fixture generation: it needs the
  # upstream CellChat package *installed* (not just checked out) and the 70 MB human-skin
  # object, neither of which a clean checkout or this job has. Its output
  # (`tests/fixtures/tutorial_netp_slice.tsv`) is committed, and both readers
  # (`gen_centrality_golden.R`, `check_centrality.R`) only ever read it behind
  # `file.exists()`, so there is nothing to regenerate and nothing that breaks by skipping.
  case "$g" in
    tests/parity/gen_tutorial_netp_slice.R) echo "== $g (one-off, skipped; output is committed)"; continue;;
  esac
  echo "== $g"
  R --vanilla -q -f "$g"
done
# The SQLite-free database exports the shim reads (`CELLCHATRS_DB`). They are *derived* from the
# same pinned `CellChatDB.<species>.rda` as the goldens, and they were the one piece `gen_fixtures.sh`
# did not produce -- so a fresh clone got past the goldens and then every gated call failed with
# "no CellChatDB export found", which reads like a broken kernel rather than a missing build step.
# Exported for every species the shim's tests and benchmarks can ask for.
for sp in human mouse; do
  echo "== export_db $sp"
  R_LIBS_USER="${R_LIBS_USER:-.rlib}" Rscript tests/parity/export_db.R "$sp" "tests/fixtures/db_$sp"
done

# f80_ref.txt is not R-generated: the oracle is a C program compiled against the platform's
# real x87 `long double` (tests/parity/gen_f80_ref.c), so the reference cannot come from R --
# it has to come from the same C compiler and the same FPU that R itself would use.
if command -v gcc >/dev/null 2>&1; then
  echo "== f80_ref"
  gcc -O2 -o /tmp/cellchatrs_f80_ref tests/parity/gen_f80_ref.c
  /tmp/cellchatrs_f80_ref > tests/fixtures/f80_ref.txt
  rm -f /tmp/cellchatrs_f80_ref
else
  echo "note: gcc not found; f80_ref.txt left as is" >&2
fi

## The standalone CLI's own fixtures. Separate from the shim's, because they are a different format
## with a different reader: the shim's goldens are read by the Rust test suite through
## `include_str!`, and the CLI's are read by a process boundary.
if [ -n "${CELLCHAT_SRC:-}" ]; then
  echo "== the CLI's input/golden pairs"
  R_LIBS_USER="${R_LIBS_USER:-.rlib}" R --vanilla -q -f tests/parity/gen_cli_fixture.R
fi

# Summary for the log: every fixture with a byte size, so a truncated regeneration shows
# up as a small file here rather than as a mysterious count mismatch in a Rust test three
# steps later. `filter_golden.txt` should hold 19 `case` records; if it holds 2, the
# generator completed yet wrote partial output, which means the R environment (not the
# generator) diverged.
echo "== fixture summary"
wc -c tests/fixtures/*.txt tests/fixtures/*.tsv 2>/dev/null | tail -25
echo "filter_golden cases: $(grep -c '^case\s' tests/fixtures/filter_golden.txt 2>/dev/null || echo 0) (expect 19)"
