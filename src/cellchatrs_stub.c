/* A translation unit so that `src/` is not empty.
 *
 * `R CMD check` warns with "Subdirectory 'src' contains no source files" when a package ships
 * `src/Makevars` and no sources. This file exists only to satisfy that check; see
 * `src/Makevars` for why the link output is renamed rather than left at the default
 * `cellchatrs.so` (which would overwrite the Rust cdylib that ./configure stages).
 *
 * It defines nothing that matters. The native entry point is `R_init_cellchatrs` in the Rust
 * crate (crates/cellchatrs/src/lib.rs).
 */

int cellchatrs_stub_unused(void);

int cellchatrs_stub_unused(void) { return 0; }
