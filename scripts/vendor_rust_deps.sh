#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/cellchatrs-vendor.XXXXXX")"
trap 'rm -rf "$tmp_dir"' EXIT HUP INT TERM

cargo vendor --locked --manifest-path "$root/src/rust/r-build/Cargo.toml" "$tmp_dir/vendor" \
  > "$tmp_dir/cargo-config.toml"
tar --sort=name --mtime='UTC 1970-01-01' --owner=0 --group=0 --numeric-owner \
  -cf - -C "$tmp_dir" vendor | xz -9e -c > "$root/src/rust/r-build/vendor.tar.xz"
(
  cd "$root/src/rust/r-build"
  sha256sum vendor.tar.xz > vendor.tar.xz.sha256
)
echo "wrote src/rust/r-build/vendor.tar.xz ($(wc -c < "$root/src/rust/r-build/vendor.tar.xz" | tr -d '[:space:]') bytes)"
