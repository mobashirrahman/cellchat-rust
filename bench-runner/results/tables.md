## Measured results

| fixture | kind | nC | K | nGenes | nLR | nboot | reps | R upstream | Rust | speedup | peak RSS | parity |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| visium | spatial | 1,073 | 8 | 433 | 134 | 100 | 5 | 17.11 s [16.78, 17.17] | 0.30 s [0.28, 0.31] | **57.22×** [54.31, 61.99] | 553 MB | divergence by design (max|Δ|=0.00177) |
| nl | scrnaseq | 2,552 | 12 | 593 | 492 | 100 | 5 | 30.27 s [29.51, 31.15] | 0.48 s [0.48, 0.51] | **62.80×** [58.21, 65.43] | 549 MB | identical |
| ls | scrnaseq | 5,011 | 12 | 593 | 691 | 100 | 5 | 44.96 s [44.56, 47.20] | 0.66 s [0.64, 0.71] | **67.92×** [62.94, 74.21] | 720 MB | identical |
| humanSkin.rda | scrnaseq | 7,563 | 12 | 1,094 | 1,583 | 100 | 5 | 101.44 s [100.82, 105.83] | 2.41 s [2.37, 2.47] | **42.14×** [40.83, 44.64] | 1270 MB | identical |
| e14 | scrnaseq | 12,179 | 13 | 582 | 832 | 100 | 5 | 83.87 s [79.79, 94.81] | 2.27 s [2.19, 2.66] | **36.93×** [29.94, 43.27] | 2062 MB | identical |
| e13 | scrnaseq | 12,951 | 11 | 582 | 835 | 100 | 5 | 78.82 s [78.02, 90.33] | 2.17 s [2.15, 2.22] | **36.32×** [35.11, 41.94] | 1880 MB | identical |
| wound.rda | scrnaseq | 21,557 | 25 | 1,101 | 1,568 | 100 | 5 | 245.34 s [229.52, 274.92] | 9.41 s [9.31, 9.63] | **26.06×** [23.83, 29.53] | 2675 MB | identical |

Protocol: 5 timed repeats after a discarded warm-up, median with a bootstrap CI, CPUs pinned with `taskset`. Raw repeats are in the JSON beside this file; the scatter vectors are dropped once their SVG is written.
