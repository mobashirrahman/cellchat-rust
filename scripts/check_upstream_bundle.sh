#!/usr/bin/env sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
BUNDLE="$ROOT/inst/upstream/CellChat-75253cd0"
PIN=75253cd0c9e68410e6e721a6d3a0419a1d7e358f

(cd "$BUNDLE" && sha256sum --check SHA256SUMS)

if [ -n "${CELLCHAT_SRC:-}" ]; then
  actual=$(git -C "$CELLCHAT_SRC" rev-parse HEAD)
  if [ "$actual" != "$PIN" ]; then
    echo "upstream checkout is $actual; expected $PIN" >&2
    exit 1
  fi
  while read -r expected path; do
    [ -n "$path" ] || continue
    cmp "$BUNDLE/$path" "$CELLCHAT_SRC/$path" || {
      echo "bundled CellChat file differs from pinned checkout: $path" >&2
      exit 1
    }
  done < "$BUNDLE/SHA256SUMS"
fi

echo "bundled CellChat files match $PIN"
