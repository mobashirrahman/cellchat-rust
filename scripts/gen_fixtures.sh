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
CC="${CELLCHAT_SRC:-/scratch/mdra00001/tmp/opencode/CellChat}"
if [ ! -d "$CC/R" ]; then
  echo "pinned upstream not found at $CC; set CELLCHAT_SRC" >&2
  exit 1
fi
for g in tests/parity/gen_*.R; do
  echo "== $g"
  R_LIBS=.rlib R --vanilla -q -f "$g"
done
# The SQLite-free database exports the shim reads (`CELLCHATRS_DB`). They are *derived* from the
# same pinned `CellChatDB.<species>.rda` as the goldens, and they were the one piece `gen_fixtures.sh`
# did not produce -- so a fresh clone got past the goldens and then every gated call failed with
# "no CellChatDB export found", which reads like a broken kernel rather than a missing build step.
# Exported for every species the shim's tests and benchmarks can ask for.
for sp in human mouse; do
  echo "== export_db $sp"
  R_LIBS=.rlib Rscript tests/parity/export_db.R "$sp" "tests/fixtures/db_$sp"
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
  R --vanilla -q -f tests/parity/gen_cli_fixture.R
fi
