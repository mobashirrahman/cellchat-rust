#!/usr/bin/env bash
# Codegen-variant parity: do different build configurations produce the same numbers?
#
# A unit test cannot observe rustc's flags, so the claim "FMA contraction is off and nothing is
# reassociated" has to be checked by building more than once and comparing. Three variants:
#
#   default       the configuration CI and the benchmarks use
#   +fma          contraction explicitly re-enabled, which must change nothing
#   native        every codegen feature this host offers, which must change nothing
#
# The check is a diff of a canonical dump of parity-critical results, not a set of pass/fail
# summaries. Comparing the numbers themselves means a variant that produced *different* wrong
# answers fails, and a variant that silently stopped running the tests fails too, which a summary
# comparison would not catch.
#
# Measured on the pinned host (8c/16t Zen3, gcc 13, rustc 1.98): all three variants produce a
# byte-identical dump. That is expected rather than lucky -- rustc does not set LLVM's `contract`
# fast-math flag on `fmul`/`fadd`, so there is no contraction to enable. The run is here so the
# property is continuously verified rather than assumed, and so a future change to the rustflags or
# to `RUSTFLAGS` in CI cannot quietly break it.
set -euo pipefail

cd "$(dirname "$0")/.."

# The parity-critical tests: the x87 oracle, the long-double special cases, the real probability
# goldens, and the FMA probes. `longdouble_special_parity` is listed separately because it is where
# a signed-zero or NaN regression from a vectorisation change would show up.
TESTS=(--test f80_vs_x87 --test no_fma --test prob_parity --test longdouble_special_parity)

OUT="${1:-target/codegen_variants}"
mkdir -p "$OUT"

# `RUSTFLAGS` replaces the config file's value rather than adding to it, so the `+fma` and `native`
# variants are set explicitly here and the default variant is set explicitly too -- otherwise
# "default" would mean "whatever the ambient environment happens to be".
run_variant() {
    local name="$1" flags="$2" file="$OUT/$1.txt"
    printf '=== %s ===\n' "$name"
    printf 'RUSTFLAGS=%s\n' "$flags"
    rm -f "$file"
    # A separate target directory per variant: sharing one would let stale objects from a previous
    # configuration satisfy the link step, and the whole point is to compare distinct builds.
    if ! CARGO_TARGET_DIR="target/codegen_$name" \
        RUSTFLAGS="$flags" \
        cargo test --release -p r-core "${TESTS[@]}" -- --test-threads=1 --nocapture \
        2>&1 | tee "$file"; then
        echo "variant $name FAILED" >&2
        return 1
    fi
    local passed
    passed="$(grep -cE '^test result: ok\. [1-9][0-9]* passed; 0 failed;' "$file" || true)"
    local expected=$(( ${#TESTS[@]} / 2 ))
    if [[ "$passed" != "$expected" ]] || grep -q '^test result: FAILED' "$file"; then
        echo "variant $name did not pass every required test binary ($passed/$expected)" >&2
        return 1
    fi
    local dump="$OUT/$name.dump"
    sed -nE 's/^.*(f80 vs x87:.*|[0-9]+ chained sums,.*|[0-9]+ cancellation chains,.*)$/\1/p' \
        "$file" | sort > "$dump"
    if [[ "$(wc -l < "$dump")" -lt 3 ]]; then
        echo "variant $name did not emit the expected parity-critical values" >&2
        return 1
    fi
}

run_variant default   "-C target-feature=-fma"
run_variant fma       "-C target-feature=+fma"
run_variant native    "-C target-cpu=native"

# Compare only the parity-critical values, not timings or target directories.
for v in fma native; do
    if ! diff -u "$OUT/default.dump" "$OUT/$v.dump" > "$OUT/default-vs-$v.diff"; then
        echo "=== $v differs from default ===" >&2
        cat "$OUT/default-vs-$v.diff" >&2
        exit 1
    fi
    echo "=== $v is bit-identical to default ==="
done

echo
echo "all three codegen variants agree; dumps in $OUT"
