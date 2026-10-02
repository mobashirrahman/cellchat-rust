#!/bin/sh
# Rebuild the Rust cdylib and reinstall the R package, **failing loudly**.
#
# The two steps must be checked separately. `cargo build | grep -c error` exits 1 when there
# are *no* errors, so `cargo ... && R CMD INSTALL ...` silently skips the install -- and R
# then loads the previous cdylib, which produces a baffling "the fix had no effect" failure
# three steps downstream. This has happened; hence the script.
set -e
cd "$(dirname "$0")/.."
echo "== cargo build --release"
if ! cargo build --release 2>&1 | tee /tmp/cellchatrs-build.log | grep -E "^(error|warning)" -A6; then
  :
fi
if grep -qE "^error" /tmp/cellchatrs-build.log; then
  echo "== cargo build FAILED; not installing" >&2
  exit 1
fi
echo "== R CMD INSTALL"
R_LIBS="$PWD/.rlib" R CMD INSTALL --no-docs --no-byte-compile --library="$PWD/.rlib" . >/tmp/cellchatrs-install.log 2>&1 || {
  echo "== R CMD INSTALL FAILED; tail of the log:" >&2
  tail -30 /tmp/cellchatrs-install.log >&2
  exit 1
}
tail -1 /tmp/cellchatrs-install.log
