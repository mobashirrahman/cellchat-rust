#!/bin/sh
# The package ships its CellChatDB exports under inst/db/ (human, mouse) so the kernel works
# out of the box, and the differential suite reads tests/fixtures/db_*. Two committed copies of
# the same export drift silently -- a human-only change to one of them would pass every gate on
# one species and fail nowhere until a user hit it. This asserts the pairs are byte-identical;
# freshness against upstream is scripts/gen_fixtures.sh's job, not this file's.
set -e
cd "$(dirname "$0")/.."
fail=0
for sp in human mouse; do
  if ! diff -r "tests/fixtures/db_$sp" "inst/db/$sp" >/dev/null 2>&1; then
    echo "MISMATCH: inst/db/$sp differs from tests/fixtures/db_$sp" >&2
    echo "  refresh with: rm -rf inst/db/$sp && cp -r tests/fixtures/db_$sp inst/db/$sp" >&2
    fail=1
  fi
done
if [ "$fail" -ne 0 ]; then exit 1; fi
echo "inst/db matches tests/fixtures for human and mouse"
