//! `r-core` — exact re-implementation of CellChat's inference kernel.
//!
//! This crate contains **no** R dependency and **no** I/O. Every public function is a
//! pure function over plain Rust data, which is what makes the parity testing in
//! `tests/` and the micro-benchmarks in `benches/` possible without an R session.
//!
//! Module map (see `docs/SEMANTICS.md` for the full contract):
//!
//! | module      | what it reproduces                                        |
//! |-------------|-----------------------------------------------------------|
//! | [`rng`]     | R's MT19937 + `sample.int` (`R_unif_index`)               |
//! | [`longdouble`] | x87 80-bit arithmetic, for R’s two-pass `mean`          |
//! | [`mathfn`]   | R’s `pnorm`, `lgammafn`, `lbeta`, `choose` (Cody/Fullerton) |
//! | [`stats`]   | `collapse::fquantile` type 7, R type-1, trimean, geometric |
//! | [`db`]      | `CellChatDB` complex / cofactor subunit tables            |
//! | [`de`]      | `computeAveExpr`, `subsetDB`, `subsetData`               |
//! | [`expr`]    | `computeExpr_LR` / `_complex` / `_coreceptor` / agonist…  |
//! | [`aggregate`] | `aggregate(matrix, list(factor), FUN)` + the `FunMean` dispatch |
//! | [`prob`]    | `computeCommunProb()` end to end — *92–95 % of runtime*  |
//! | [`net`]     | `aggregateNet`, `subsetCommunication`, the melt order   |
//! | [`wilcox`]   | `rank`, `wilcox.test`, `p.adjust`, `mean.fxn`              |

/// The crate version, so a driver built against this source can report which numerics it has
/// rather than only its own version. A CLI that prints one number and calls it "the version" is
/// ambiguous exactly when it matters.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub mod aggregate;
pub mod centrality;
pub mod db;
pub mod de;
pub mod expr;
pub mod filter;
pub mod longdouble;
pub mod mathfn;
pub mod net;
pub mod pathway;
pub mod prob;
pub mod ranknet;
pub mod rng;
pub mod snn;
pub mod spatial;
pub mod stats;
pub mod subset;
pub mod wilcox;

/// How closely results must match upstream R.
///
/// Recorded in the CI-published `parity.json` so the claim is machine-checkable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParityRung {
    /// byte-for-byte identical
    Exact,
    /// within 1 ulp of upstream
    OneUlp,
    /// within a stated relative tolerance
    Approx,
}
