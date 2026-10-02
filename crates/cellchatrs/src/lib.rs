//! extendr bindings. **Argument marshalling only** — no numerics live here, so the
//! binding layer cannot silently diverge from `r-core`.
//!
//! ## Wiring notes (all three verified empirically on this machine)
//!
//! 1. The umbrella crate `extendr` is not resolvable from the crates.io sparse index
//!    (HTTP 404); depend on `extendr-api` + `extendr-macros` directly.
//! 2. `extendr_module!` changed syntax in 0.9 — it now takes a braced list:
//!    `extendr_module! { mod cellchatrs; fn compute_commun_prob; }`.
//!    Note the module name must match the crate name, because R derives the init
//!    symbol name from the DLL name.
//! 3. The generated init symbol is `R_init_<mod>_extendr`, which R's `dyn.load()` does
//!    *not* auto-invoke. We therefore register routines from the metadata table in an
//!    explicit `init` entry point rather than relying on autoload; see `R/zzz.R`.

use extendr_api::prelude::*;

// The `#[extendr]`-annotated functions live at the crate root, not in a submodule, because
// the macro emits a **private** `meta__<name>` next to each annotated function and
// `extendr_module!` expands to `crate::meta__<name>` paths. In a submodule that path is
// inaccessible from the crate root; in the crate root it is in scope. `include!` rather than
// `mod`, so the items land in the right module.
include!("prob.rs");
include!("wilcox.rs");

/// Requested rayon thread count, read once when the global pool is first built.
///
/// `rayon::ThreadPoolBuilder::build_global()` can only be called **once per process**,
/// so a plain `set_num_threads()` that rebuilds the global pool silently no-ops if the
/// kernel has already run. Instead we record the request in an atomic and build the
/// global pool lazily on first use, so "configure then run" always takes effect.
/// (This was a real bug found while validating the wiring: the first version built the
/// pool in `.onLoad` and every later `set_num_threads()` was a no-op.)
static REQUESTED_THREADS: std::sync::atomic::AtomicUsize =
    std::sync::atomic::AtomicUsize::new(usize::MAX);

/// Build the global rayon pool on first use, honouring any prior `set_num_threads`.
fn ensure_pool() {
    use std::sync::Once;
    static INIT: Once = Once::new();
    INIT.call_once(|| {
        let requested = REQUESTED_THREADS.load(std::sync::atomic::Ordering::Relaxed);
        let n = if requested == usize::MAX {
            std::thread::available_parallelism()
                .map(|v| v.get())
                .unwrap_or(1)
        } else {
            requested
        };
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global();
    });
}

/// Record a thread-count request **without** building the pool. Safe to call at any time;
/// takes effect if the pool has not been built yet. This is what `.onLoad` uses.
#[extendr]
fn request_num_threads(n: usize) -> usize {
    REQUESTED_THREADS.store(n, std::sync::atomic::Ordering::Relaxed);
    n
}

/// Size the rayon pool. Only effective before the first kernel call; returns the size
/// actually in force afterwards.
#[extendr]
fn set_num_threads(n: usize) -> usize {
    request_num_threads(n);
    ensure_pool();
    rayon::current_num_threads()
}

#[extendr]
fn get_num_threads() -> usize {
    ensure_pool();
    rayon::current_num_threads()
}

extendr_module! {
    mod cellchatrs;
    fn request_num_threads;
    fn set_num_threads;
    fn get_num_threads;
    fn db_load;
    fn average_expression;
    fn compute_commun_prob;
    fn resolve_cofactor_rows;
    fn match_mean_type;
    fn aggregate_net;
    fn subset_communication;
    fn subset_communication_deg;
    fn ranknet_information_flow;
    fn compute_ave_expr;
    fn subset_data_gene_use;
    fn subset_db_by_annotation;
    fn identify_over_expressed_genes;
    fn identify_over_expressed_genes_dataset;
    fn expressed_in_cells;
    fn compute_commun_prob_pathway;
    fn aggregate_net_filtered;
    fn filter_communication;
    fn spatial_trimmed_mean;
    fn spatial_fdist;
    fn spatial_region_distance;
    fn ranknet_comparison_flow;
    fn ranknet_pairwise_orders;
    fn centrality_deterministic;
}

/// R's `dyn.load()` auto-invokes `R_init_<dll-filename>` — note **no** `_extendr` suffix.
/// extendr 0.9 emits `R_init_<mod>_extendr`, so without this shim nothing is ever
/// registered and R reports "could not find function". Verified empirically.
#[no_mangle]
#[allow(non_snake_case)]
pub extern "C" fn R_init_cellchatrs(info: *mut extendr_api::DllInfo) {
    R_init_cellchatrs_extendr(info);
}
