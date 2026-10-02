//! `CellChatDB`: the ligand–receptor database, parsed into flat index-addressed tables.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md`). The database is not numerically
//! subtle, but three of its *orderings* are load-bearing and are easy to get subtly
//! wrong, so each is called out below and pinned by a test.
//!
//! ## Why a text export rather than reading the `.rda`
//!
//! `r-core` must have no R dependency and no I/O so its numerics can be differentially
//! tested without an R session. `tests/parity/export_db.R` flattens the pinned
//! `CellChatDB.{species}.rda` into four TSV files plus a manifest, and the manifest's
//! MD5s let [`Database::verify_manifest`] refuse a stale fixture.
//!
//! ## The three load-bearing orderings
//!
//! 1. **Interaction order is the annotation factor order.** `subsetData`
//!    (`R/utilities.R:319`) does
//!    `interaction_input[order(interaction_input$annotation), ]` while `annotation` is a
//!    *factor* with levels `Secreted Signaling, ECM-Receptor, Non-protein Signaling,
//!    Cell-Cell Contact`, then converts back to character. `computeCommunProb` derives
//!    `nLR1` from the boundary between the first three and the last, so this ordering
//!    determines which L-R pairs are evaluated in the contact-dependent regime. The
//!    export applies the sort, and verifies that R's `order()` is a *stable* partition
//!    (it is, on R 4.3.3) so the within-level order is preserved.
//! 2. **Subunit order within a complex is `subunit_1 .. subunit_5`,** in that order, with
//!    empty cells dropped. `computeExpr_complex` takes a *geometric* mean, which is
//!    order-independent — but `unlist()` on the data frame is column-major, so a naive
//!    port that reorders subunits changes nothing numerically while breaking any future
//!    hash/equality check. Kept faithful anyway.
//! 3. **`extractGeneSubset` output order** is `intersect(geneSet, symbols)` followed by
//!    the complex subunits, deduplicated. It does not reach the kernel, because
//!    `subsetData` selects rows with `rownames(data) %in% gene.use` and therefore
//!    inherits the *expression matrix's* row order. It is preserved because
//!    `identifyOverExpressedGenes` and `computeAveExpr(features = ...)` do depend on it.

use std::collections::HashMap;
use std::fmt;
use std::path::Path;

/// The annotation factor levels, in the order `subsetData` uses.
///
/// The first three are diffusion-mediated; `Cell-Cell Contact` is contact-dependent.
/// [`Database::nlr1`] is the last index belonging to the first three (0 if none).
pub const ANNOTATION_LEVELS: [&str; 4] = [
    "Secreted Signaling",
    "ECM-Receptor",
    "Non-protein Signaling",
    "Cell-Cell Contact",
];

/// A ligand–receptor interaction, with the optional regulators CellChatDB attaches.
///
/// Empty strings are meaningful throughout: they are how R encodes `NA` in the `.rda`
/// and how upstream tests for "no agonist" (`!is.na(x) & x != ""`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Interaction {
    pub name: String,
    pub pathway: String,
    pub ligand: String,
    pub receptor: String,
    pub agonist: String,
    pub antagonist: String,
    pub co_activator: String,
    pub co_inhibitor: String,
    pub annotation: String,
}

impl Interaction {
    /// True when the pair has a non-empty, present agonist.
    pub fn has_agonist(&self) -> bool {
        !self.agonist.is_empty()
    }

    /// True when the pair has a non-empty antagonist.
    pub fn has_antagonist(&self) -> bool {
        !self.antagonist.is_empty()
    }

    /// True when the pair is contact-dependent (i.e. the last annotation level).
    pub fn is_contact_dependent(&self) -> bool {
        self.annotation == ANNOTATION_LEVELS[3]
    }

    /// True when the pair is diffusion-mediated (one of the first three levels).
    pub fn is_diffusion_mediated(&self) -> bool {
        self.annotation == ANNOTATION_LEVELS[0]
            || self.annotation == ANNOTATION_LEVELS[1]
            || self.annotation == ANNOTATION_LEVELS[2]
    }
}

/// Parse errors are reported rather than silently ignored: a mis-parsed database would
/// produce a plausible but wrong network, which is the failure mode this whole project
/// exists to prevent.
#[derive(Debug)]
pub enum DbError {
    Io(std::io::Error),
    Malformed {
        file: &'static str,
        line: usize,
        msg: String,
    },
    UnknownAnnotation(String),
    Stale {
        what: &'static str,
        expected: String,
        found: String,
    },
}

impl fmt::Display for DbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DbError::Io(e) => write!(f, "cannot read CellChatDB export: {e}"),
            DbError::Malformed { file, line, msg } => {
                write!(f, "{file}:{line}: {msg}")
            }
            DbError::UnknownAnnotation(a) => write!(
                f,
                "annotation {a:?} is not one of {:?}; it would change the \
                 contact/diffusion split and therefore nLR1",
                ANNOTATION_LEVELS
            ),
            DbError::Stale {
                what,
                expected,
                found,
            } => write!(
                f,
                "CellChatDB export is stale: {what} has MD5 {found}, manifest expects \
                 {expected}. Re-run tests/parity/export_db.R against the pinned commit."
            ),
        }
    }
}

impl std::error::Error for DbError {}

impl From<std::io::Error> for DbError {
    fn from(e: std::io::Error) -> Self {
        DbError::Io(e)
    }
}

type Result<T> = std::result::Result<T, DbError>;

/// The full database: interactions, the two subunit tables, and the official-symbol
/// oracle, plus the name→index maps that turn string lookups in the inner loop into
/// array indexing.
#[derive(Clone, Debug, Default)]
pub struct Database {
    pub species: String,
    pub interactions: Vec<Interaction>,
    /// Complex name → its `subunit_*` cells, **empties preserved**.
    pub complexes: HashMap<String, Vec<String>>,
    /// Cofactor name → its `cofactor*` cells, **empties preserved**.
    pub cofactors: HashMap<String, Vec<String>>,
    /// Number of `subunit_*` / `cofactor*` columns in the source tables. Needed because
    /// the empty cells carry the column positions (see [`Database::extract_gene_subset`]).
    pub n_subunit_cols: usize,
    pub n_cofactor_cols: usize,
    /// Official gene symbols from `geneInfo$Symbol`.
    pub symbols: HashMap<String, ()>,
    /// Complex names in **source row order** (`complex$complex_name`).
    ///
    /// `complexes` is a `HashMap`, so this is the only ordered view. Anything that must
    /// reproduce R's `unlist(complex_input[, subunit_cols])` -- which walks rows in table
    /// order -- has to iterate this instead. Getting this wrong is invisible in a
    /// per-pair test and loudly wrong in a fixture that lays out an expression matrix by
    /// "every subunit in the database", which is exactly what the `computeExpr` corpus
    /// does: the gene *set* is unchanged, only the order, so every name still resolves
    /// and every value is still plausible, but the column-major draw lands on different
    /// genes.
    pub complex_names: Vec<String>,
    /// Cofactor names in **source row order** (`cofactor$cofactor_name`).
    pub cofactor_names: Vec<String>,
    /// Complex/cofactor name → row index, for O(1) lookup in the hot loop.
    complex_index: HashMap<String, usize>,
    cofactor_index: HashMap<String, usize>,
    symbol_index: HashMap<String, usize>,
    /// 0-based index of the last diffusion-mediated interaction, i.e. upstream's
    /// 1-based `nLR1`. `0` means "all contact-dependent".
    pub nlr1: usize,
}

impl Database {
    /// Build a database from parts, filling in the private name -> index maps.
    ///
    /// Needed whenever the complex/cofactor tables and the symbol list come from somewhere
    /// other than a pinned export directory -- a `CellChat` object's own `@DB` slot, or a
    /// corpus fixture. `complex_names` must be in **source row order**, because
    /// `extract_gene_subset` reproduces R's `unlist(complex_input[, subunit_cols])`, which
    /// walks rows in table order; the `complexes` map alone cannot express that.
    #[allow(clippy::too_many_arguments)]
    pub fn from_parts(
        species: impl Into<String>,
        interactions: Vec<Interaction>,
        complexes: HashMap<String, Vec<String>>,
        cofactors: HashMap<String, Vec<String>>,
        complex_names: Vec<String>,
        cofactor_names: Vec<String>,
        n_subunit_cols: usize,
        n_cofactor_cols: usize,
        symbols: HashMap<String, ()>,
    ) -> Self {
        let mut db = Database {
            species: species.into(),
            interactions,
            complexes,
            cofactors,
            complex_names,
            cofactor_names,
            n_subunit_cols,
            n_cofactor_cols,
            symbols,
            complex_index: HashMap::new(),
            cofactor_index: HashMap::new(),
            symbol_index: HashMap::new(),
            // Not derivable from the parts; only `computeCommunProb`'s spatial tail reads
            // it, and a database built here is never used for that. `0` means "no split",
            // which is the RNA behaviour.
            nlr1: 0,
        };
        db.reindex();
        db
    }

    /// Recompute the private lookup tables from the public ones. Called by [`Self::load`] and
    /// [`Self::from_parts`]; a caller that mutates the public fields must call it again.
    pub fn reindex(&mut self) {
        self.symbol_index = self
            .symbols
            .keys()
            .enumerate()
            .map(|(i, s)| (s.clone(), i))
            .collect();
        self.complex_index = self
            .complex_names
            .iter()
            .enumerate()
            .map(|(i, k)| (k.clone(), i))
            .collect();
        self.cofactor_index = self
            .cofactor_names
            .iter()
            .enumerate()
            .map(|(i, k)| (k.clone(), i))
            .collect();
    }

    /// Parse an export directory produced by `tests/parity/export_db.R`.
    pub fn load(dir: &Path) -> Result<Self> {
        let manifest = read_manifest(dir)?;
        let interactions_text = read_to_string(dir, "interactions.tsv")?;
        let complexes_text = read_to_string(dir, "complexes.tsv")?;
        let cofactors_text = read_to_string(dir, "cofactors.tsv")?;
        let symbols_text = read_to_string(dir, "symbols.txt")?;

        verify_hash(
            "interactions.tsv",
            &interactions_text,
            &manifest,
            "md5_interactions",
        )?;
        verify_hash("complexes.tsv", &complexes_text, &manifest, "md5_complexes")?;
        verify_hash("cofactors.tsv", &cofactors_text, &manifest, "md5_cofactors")?;
        verify_hash("symbols.txt", &symbols_text, &manifest, "md5_symbols")?;

        let symbols: HashMap<String, ()> = symbols_text
            .lines()
            .filter(|l| !l.is_empty())
            .map(|l| (l.to_string(), ()))
            .collect();

        let (complex_names, complexes) = parse_subunit_table(&complexes_text, "complexes.tsv")?;
        let (cofactor_names, cofactors) = parse_subunit_table(&cofactors_text, "cofactors.tsv")?;
        let n_subunit_cols = complexes.values().map(|v| v.len()).max().unwrap_or(0);
        let n_cofactor_cols = cofactors.values().map(|v| v.len()).max().unwrap_or(0);
        let mut db = Database {
            species: manifest.get("species").cloned().unwrap_or_default(),
            interactions: parse_interactions(&interactions_text)?,
            complexes,
            cofactors,
            complex_names,
            cofactor_names,
            n_subunit_cols,
            n_cofactor_cols,
            symbols,
            ..Default::default()
        };

        // Dense name -> row indices for the hot loop. `HashMap<String, usize>` costs a
        // string hash per lookup; these let `expr.rs` resolve a subunit list once per
        // L-R pair and then work entirely in `u32` indices.
        db.reindex();

        db.nlr1 = db
            .interactions
            .iter()
            .rposition(|i| i.is_diffusion_mediated())
            .map(|i| i + 1)
            .unwrap_or(0);
        if let Some(expected) = manifest.get("nlr1").and_then(|v| v.parse::<usize>().ok()) {
            if expected != db.nlr1 {
                return Err(DbError::Stale {
                    what: "nlr1",
                    expected: expected.to_string(),
                    found: db.nlr1.to_string(),
                });
            }
        }
        Ok(db)
    }

    /// Check the manifest against the pinned upstream commit.
    pub fn verify_manifest(&self, dir: &Path) -> Result<()> {
        let m = read_manifest(dir)?;
        const PIN: &str = "75253cd0c9e68410e6e721a6d3a0419a1d7e358f";
        if let Some(c) = m.get("upstream_commit") {
            if c != PIN {
                return Err(DbError::Stale {
                    what: "upstream_commit",
                    expected: PIN.to_string(),
                    found: c.to_string(),
                });
            }
        }
        Ok(())
    }

    /// Number of interactions, i.e. upstream's `nLR`.
    pub fn n_lr(&self) -> usize {
        self.interactions.len()
    }

    /// Is `name` an official gene symbol?
    ///
    /// The `.rda` contains exactly one `NA` in `geneInfo$Symbol`; NA is not a symbol,
    /// because `%in%` only matches an `NA` needle against an `NA` table entry and needles
    /// are never NA in this code path.
    pub fn is_symbol(&self, name: &str) -> bool {
        self.symbols.contains_key(name)
    }

    /// 0-based index of a symbol, for O(1) gene lookups in the hot loop.
    pub fn symbol_id(&self, name: &str) -> Option<usize> {
        self.symbol_index.get(name).copied()
    }

    /// 0-based index of a complex row, or `None` if the name is not in the table.
    pub fn complex_id(&self, name: &str) -> Option<usize> {
        self.complex_index.get(name).copied()
    }

    /// 0-based index of a cofactor row, or `None` if the name is not in the table.
    pub fn cofactor_id(&self, name: &str) -> Option<usize> {
        self.cofactor_index.get(name).copied()
    }

    /// Every symbol, in the order they were loaded.
    pub fn symbol_list(&self) -> Vec<&str> {
        let mut v: Vec<&str> = self.symbols.keys().map(|s| s.as_str()).collect();
        v.sort_unstable();
        v
    }

    /// Raw `subunit_*` cells for a complex, empties included. Order is `subunit_1..`.
    pub fn complex_cells(&self, name: &str) -> Option<&[String]> {
        self.complexes.get(name).map(|v| v.as_slice())
    }

    /// Raw `cofactor*` cells for a cofactor, empties included.
    pub fn cofactor_cells(&self, name: &str) -> Option<&[String]> {
        self.cofactors.get(name).map(|v| v.as_slice())
    }

    /// Non-empty subunits of a complex, in `subunit_*` order. Convenient for reading;
    /// the kernel uses the raw cells so it can reproduce `unlist`'s column-major order.
    pub fn complex_subunits(&self, name: &str) -> Option<Vec<&str>> {
        self.complexes.get(name).map(|v| {
            v.iter()
                .filter(|s| !s.is_empty())
                .map(|s| s.as_str())
                .collect()
        })
    }

    /// Non-empty subunits of a cofactor, in `cofactor*` order.
    pub fn cofactor_subunits(&self, name: &str) -> Option<Vec<&str>> {
        self.cofactors.get(name).map(|v| {
            v.iter()
                .filter(|s| !s.is_empty())
                .map(|s| s.as_str())
                .collect()
        })
    }

    /// Port of `extractGeneSubset` (`R/database.R:221`).
    ///
    /// ```r
    /// complex <- geneSet[which(geneSet %in% geneIfo$Symbol == "FALSE")]
    /// geneSet <- intersect(geneSet, geneIfo$Symbol)
    /// complexsubunits <- select(complex_input[match(complex, rownames(complex_input), nomatch=0),],
    ///                          starts_with("subunit"))
    /// complex <- intersect(complex, rownames(complexsubunits))
    /// complexsubunitsV <- unique(unlist(complexsubunits)[unlist(complexsubunits) != ""])
    /// unique(c(geneSet, complexsubunitsV))
    /// ```
    ///
    /// The `== "FALSE"` comparison is the load-bearing oddity. `geneSet %in% geneIfo$Symbol`
    /// is a logical vector, and R coerces it to character before comparing, so
    /// `as.character(FALSE) == "FALSE"` is TRUE. The branch therefore selects the
    /// members of `geneSet` that are **absent** from `geneInfo$Symbol` — i.e. the complex
    /// names, which are not official symbols. The shipped DB has zero literal `"FALSE"`
    /// symbols, so reading it literally (as "the names that equal FALSE") would silently
    /// drop every complex and every one of its subunits.
    ///
    /// Order of the result is `intersect(...)` order followed by the subunits, which is
    /// what upstream returns and what `identifyOverExpressedGenes` observes.
    pub fn extract_gene_subset(&self, gene_set: &[String]) -> Vec<String> {
        // `complex`: names in `gene_set` that are NOT official symbols.
        let complex: Vec<&String> = gene_set.iter().filter(|g| !self.is_symbol(g)).collect();
        // `geneSet <- intersect(geneSet, geneIfo$Symbol)`: R's `intersect` preserves the
        // order of its *first* argument, so this is the input order filtered to symbols.
        let mut out: Vec<String> = Vec::with_capacity(gene_set.len());
        let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for g in gene_set {
            if self.is_symbol(g) && seen.insert(g.as_str()) {
                out.push(g.clone());
            }
        }
        // `complex <- intersect(complex, rownames(complexsubunits))`: only names with a
        // row in the complex table contribute, and `match()` keeps the input order.
        let rows: Vec<&Vec<String>> = complex
            .iter()
            .filter_map(|n| self.complexes.get(n.as_str()))
            .collect();
        // `complexsubunitsV <- unique(unlist(complexsubunits))`, and `unlist` on a
        // multi-row data frame is **COLUMN-major**: every row's `subunit_1`, then every
        // row's `subunit_2`, and so on. Emitting row-major instead reorders the subunits.
        // Real example from the human DB: `IL12AB` is `IL12A, <empty>, IL12B`, so
        // dropping the empty cell would move IL12B into the `subunit_2` slot and shift
        // every later subunit.
        for col in 0..self.n_subunit_cols {
            for r in &rows {
                if let Some(v) = r.get(col) {
                    if !v.is_empty() && seen.insert(v.as_str()) {
                        out.push(v.clone());
                    }
                }
            }
        }
        out
    }

    /// Port of `extractGene` (`R/database.R:174`): the full gene set for a database,
    /// i.e. ligands and receptors expanded through complexes, plus all cofactor subunits.
    pub fn extract_gene(&self) -> Vec<String> {
        let ligs: Vec<String> = self.distinct(|i| &i.ligand);
        let recs: Vec<String> = self.distinct(|i| &i.receptor);
        let mut gene_lr = self.extract_gene_subset(&ligs);
        gene_lr.extend(self.extract_gene_subset(&recs));

        // `cofactor <- unique(c(agonist, antagonist, co_A_receptor, co_I_receptor))`.
        //
        // NB the concatenation is **column-wise**: R's `c()` over four vectors appends
        // every interaction's agonist, then every antagonist, then every co-activator,
        // then every co-inhibitor. Iterating the rows and reading the four fields per row
        // interleaves them instead, which changes the order of the subunits that follow
        // (the cofactor list is deduplicated but not sorted) and so changes the output
        // order of `extractGene`.
        let mut cof: Vec<&str> = Vec::new();
        const FIELDS: [u8; 4] = [0, 1, 2, 3]; // agonist, antagonist, co_A, co_I
        for field in FIELDS {
            for i in &self.interactions {
                let v: &String = match field {
                    0 => &i.agonist,
                    1 => &i.antagonist,
                    2 => &i.co_activator,
                    _ => &i.co_inhibitor,
                };
                // `cofactor <- cofactor[cofactor != ""]` happens after `unique`, but
                // dropping empties during the scan is equivalent: they are never
                // referenced again and contribute no subunits.
                if !v.is_empty() && !cof.contains(&v.as_str()) {
                    cof.push(v.as_str());
                }
            }
        }
        // Owned keys: `gene_lr` is a local, so borrowing from it into the `seen` set
        // would not outlive the loop.
        let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        let mut out: Vec<String> = Vec::with_capacity(gene_lr.len());
        for g in gene_lr {
            if seen.insert(g.clone()) {
                out.push(g);
            }
        }
        // Column-major, for the same `unlist` reason as `extract_gene_subset`.
        let rows: Vec<&Vec<String>> = cof.iter().filter_map(|n| self.cofactors.get(*n)).collect();
        for col in 0..self.n_cofactor_cols {
            for r in &rows {
                if let Some(v) = r.get(col) {
                    if !v.is_empty() && seen.insert(v.to_string()) {
                        out.push(v.clone());
                    }
                }
            }
        }
        out
    }

    fn distinct<F: Fn(&Interaction) -> &String>(&self, f: F) -> Vec<String> {
        let mut seen = std::collections::HashSet::new();
        let mut out = Vec::new();
        for i in &self.interactions {
            let v = f(i);
            if seen.insert(v.as_str()) {
                out.push(v.clone());
            }
        }
        out
    }

    /// Interactions in the diffusion-mediated regime, in export order.
    pub fn diffusion_mediated(&self) -> impl Iterator<Item = &Interaction> {
        self.interactions
            .iter()
            .filter(|i| i.is_diffusion_mediated())
    }

    /// Interactions in the contact-dependent regime, in export order.
    pub fn contact_dependent(&self) -> impl Iterator<Item = &Interaction> {
        self.interactions
            .iter()
            .filter(|i| i.is_contact_dependent())
    }
}

// ---------------------------------------------------------------------- parsing

fn read_to_string(dir: &Path, name: &'static str) -> Result<String> {
    std::fs::read_to_string(dir.join(name)).map_err(|e| DbError::Malformed {
        file: name,
        line: 0,
        msg: format!("{e}"),
    })
}

type Manifest = HashMap<String, String>;

fn read_manifest(dir: &Path) -> Result<Manifest> {
    let text = read_to_string(dir, "MANIFEST.tsv")?;
    let mut m = HashMap::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
    {
        if let Some((k, v)) = line.split_once('\t') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    Ok(m)
}

/// MD5, matching `tools::md5sum`.
///
/// The export manifest is a *staleness* check, not an integrity guarantee, so a
/// self-contained MD5 is proportionate; pulling in `md-5`/`digest` for a stronger hash
/// would add a dependency to satisfy nothing.
fn md5_hex(bytes: &[u8]) -> String {
    const S: [u32; 64] = [
        7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 5, 9, 14, 20, 5, 9, 14, 20, 5,
        9, 14, 20, 5, 9, 14, 20, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 6, 10,
        15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21,
    ];
    let mut k = [0u32; 64];
    for (i, slot) in k.iter_mut().enumerate() {
        *slot = ((i as f64 + 1.0).sin().abs() * 4294967296.0) as u32;
    }
    let (mut a0, mut b0, mut c0, mut d0) =
        (0x67452301u32, 0xefcdab89u32, 0x98badcfeu32, 0x10325476u32);
    let mut msg = bytes.to_vec();
    let bitlen = (bytes.len() as u64).wrapping_mul(8);
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bitlen.to_le_bytes());
    for chunk in msg.chunks(64) {
        let mut m = [0u32; 16];
        for (i, slot) in m.iter_mut().enumerate() {
            *slot = u32::from_le_bytes(chunk[4 * i..4 * i + 4].try_into().unwrap());
        }
        let (mut a, mut b, mut c, mut d) = (a0, b0, c0, d0);
        for i in 0..64 {
            let (f, g) = match i / 16 {
                0 => ((b & c) | (!b & d), i),
                1 => ((d & b) | (!d & c), (5 * i + 1) % 16),
                2 => (b ^ c ^ d, (3 * i + 5) % 16),
                _ => (c ^ (b | !d), (7 * i) % 16),
            };
            let f2 = f.wrapping_add(a).wrapping_add(k[i]).wrapping_add(m[g]);
            a = d;
            d = c;
            c = b;
            b = b.wrapping_add(f2.rotate_left(S[i]));
        }
        a0 = a0.wrapping_add(a);
        b0 = b0.wrapping_add(b);
        c0 = c0.wrapping_add(c);
        d0 = d0.wrapping_add(d);
    }
    let mut out = String::with_capacity(32);
    for v in [a0, b0, c0, d0] {
        for byte in v.to_le_bytes() {
            out.push_str(&format!("{byte:02x}"));
        }
    }
    out
}

fn verify_hash(name: &'static str, text: &str, manifest: &Manifest, key: &str) -> Result<()> {
    // `writeLines` appends a trailing newline; the manifest hashed the file on disk, so
    // compare against that rather than against a re-serialised copy.
    if let Some(expected) = manifest.get(key) {
        let found = md5_hex(text.as_bytes());
        if &found != expected {
            return Err(DbError::Stale {
                what: name,
                expected: expected.clone(),
                found,
            });
        }
    }
    Ok(())
}

fn parse_interactions(text: &str) -> Result<Vec<Interaction>> {
    let mut out = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != 9 {
            return Err(DbError::Malformed {
                file: "interactions.tsv",
                line: i + 1,
                msg: format!("expected 9 tab-separated fields, got {}", f.len()),
            });
        }
        let annotation = f[8].to_string();
        if !ANNOTATION_LEVELS.contains(&annotation.as_str()) {
            return Err(DbError::UnknownAnnotation(annotation));
        }
        out.push(Interaction {
            name: f[0].to_string(),
            pathway: f[1].to_string(),
            ligand: f[2].to_string(),
            receptor: f[3].to_string(),
            agonist: f[4].to_string(),
            antagonist: f[5].to_string(),
            co_activator: f[6].to_string(),
            co_inhibitor: f[7].to_string(),
            annotation,
        });
    }
    Ok(out)
}

/// Parse a `name<TAB>cell<TAB>cell...` table, **keeping empty cells**.
///
/// The empties are load-bearing: they carry the column positions that `unlist`'s
/// column-major flattening depends on (see [`Database::extract_gene_subset`]). Trailing
/// tabs survive because `str::split('\t')` does not drop trailing empty fields.
/// A row with no cells at all is legal and yields an empty list rather than being
/// dropped, so a lookup for that name still succeeds.
fn parse_subunit_table(
    text: &str,
    file: &'static str,
) -> Result<(Vec<String>, HashMap<String, Vec<String>>)> {
    let mut order = Vec::new();
    let mut out = HashMap::new();
    for (i, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.is_empty() || f[0].is_empty() {
            return Err(DbError::Malformed {
                file,
                line: i + 1,
                msg: "row has no name field".into(),
            });
        }
        out.insert(
            f[0].to_string(),
            f[1..].iter().map(|s| s.to_string()).collect(),
        );
        order.push(f[0].to_string());
    }
    Ok((order, out))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn md5_matches_known_vectors() {
        // RFC 1321 test suite.
        for (input, want) in [
            ("", "d41d8cd98f00b204e9800998ecf8427e"),
            ("a", "0cc175b9c0f1b6a831c399e269772661"),
            ("abc", "900150983cd24fb0d6963f7d28e17f72"),
            ("message digest", "f96b697d7cb7938d525a2f31aaf161d0"),
            (
                "abcdefghijklmnopqrstuvwxyz",
                "c3fcd3d76192e4007dfb496cca67e13b",
            ),
            (
                "12345678901234567890123456789012345678901234567890123456789012345678901234567890",
                "57edf4a22be3c955ac49da2e2107b67a",
            ),
        ] {
            assert_eq!(md5_hex(input.as_bytes()), want, "md5({input:?})");
        }
    }

    #[test]
    fn interaction_flags_match_the_annotations() {
        let i = Interaction {
            name: "X".into(),
            pathway: "P".into(),
            ligand: "A".into(),
            receptor: "B".into(),
            agonist: String::new(),
            antagonist: "ant".into(),
            co_activator: String::new(),
            co_inhibitor: String::new(),
            annotation: ANNOTATION_LEVELS[3].into(),
        };
        assert!(!i.has_agonist() && i.has_antagonist());
        assert!(i.is_contact_dependent() && !i.is_diffusion_mediated());
    }

    #[test]
    fn empty_and_nan_regulators_are_both_absent() {
        let mk = |a: &str, an: &str| Interaction {
            name: "X".into(),
            pathway: "P".into(),
            ligand: "A".into(),
            receptor: "B".into(),
            agonist: a.into(),
            antagonist: an.into(),
            co_activator: String::new(),
            co_inhibitor: String::new(),
            annotation: ANNOTATION_LEVELS[0].into(),
        };
        assert!(!mk("", "").has_agonist());
        assert!(!mk("", "").has_antagonist());
        assert!(mk("ag", "").has_agonist());
    }

    #[test]
    fn unknown_annotation_is_rejected() {
        let line = "N\tP\tL\tR\t\t\t\t\tBogus Signaling\n";
        assert!(matches!(
            parse_interactions(line),
            Err(DbError::UnknownAnnotation(a)) if a == "Bogus Signaling"
        ));
    }

    #[test]
    fn short_row_is_rejected() {
        let line = "N\tP\tL\n";
        assert!(matches!(
            parse_interactions(line),
            Err(DbError::Malformed { line: 1, .. })
        ));
    }
}
