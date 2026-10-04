/* A translation unit so that `src/` is not empty.
 *
 * `R CMD check` warns with "Subdirectory 'src' contains no source files" when a package ships
 * `src/Makevars` and no sources. This file exists only to satisfy that check; see
 * `src/Makevars` for why the link output is named `CellChat.so` (which is what `$(SHLIB)`
 * resolves to for this package) rather than left to overwrite the Rust cdylib that
 * ./configure stages.
 *
 * It defines nothing that matters. The native entry point is `R_init_CellChat` in the Rust
 * crate (src/rust/crates/cellchatrs/src/lib.rs).
 */

int CellChat_stub_unused(void);

int CellChat_stub_unused(void) { return 0; }
