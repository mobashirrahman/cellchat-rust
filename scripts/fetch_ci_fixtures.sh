#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
mkdir -p data

fetch() {
  local name="$1" file_id="$2" bytes="$3" md5="$4" path="data/$1"
  curl -fLsS --retry 3 --retry-delay 2 \
    "https://ndownloader.figshare.com/files/$file_id" -o "$path"
  local actual_bytes
  actual_bytes="$(wc -c < "$path" | tr -d '[:space:]')"
  if [[ "$actual_bytes" != "$bytes" ]]; then
    echo "fixture $name has $actual_bytes bytes; expected $bytes" >&2
    exit 1
  fi
  if ! printf '%s  %s\n' "$md5" "$path" | md5sum --check --status; then
    echo "fixture $name failed its Figshare MD5 check" >&2
    exit 1
  fi
  echo "fetched $name ($bytes bytes, MD5 verified)"
}

case "${1:-}" in
  tutorial)
    fetch humanSkin.rda 42997198 13426240 9d31f4c88be64912b937b5d6a2cb6c2d
    ;;
  parity)
    fetch humanSkin.rda 42997198 13426240 9d31f4c88be64912b937b5d6a2cb6c2d
    fetch wound.rda 38838357 57217634 e2bf9042fe1fb41894c55ec5d57475d6
    fetch visium.rds 43061620 2246226 4caf7f60423f75ab7870fb7a30c6c412
    ;;
  *)
    echo "usage: $0 tutorial|parity" >&2
    exit 2
    ;;
esac
