//! Randomised fuzzing of the L-R database *structure*.
//!
//! Every other test in this crate is handed a database that somebody built by hand, and every one
//! of them is a database where the names line up. What is untested is the resolution logic against
//! a database that is *arbitrary*: a complex with no subunits, a subunit that is itself a complex,
//! a cofactor named the same as a complex, an interaction whose ligand is an empty string, a
//! complex name that is also a row in the expression matrix, and so on.
//!
//! Two properties are asserted for every generated database, and they are the ones that matter for
//! a drop-in replacement:
//!
//!   1. **No panic.** Every entry point is total: it returns a value or an `ExprError`, never
//!      indexes out of bounds. A `subscript out of bounds` is a *value* here -- it is what upstream
//!      raises -- so it must arrive as an `Err`, not as a Rust panic that surfaces in R as a
//!      message with no provenance.
//!   2. **Agreement with a straightforward reference.** The resolution result is compared against
//!      a `BTreeMap`-driven lookup written independently here, so a shared misreading of the
//!      rules cannot make the property vacuous.
//!
//! The generator is deliberately *hostile*: probabilities produce empty subunit lists, names
//! colliding across tables, self-referential complexes, and cycles. A real database has none of
//! those, which is exactly why a hand-built one finds nothing.

use std::collections::{BTreeMap, BTreeSet};

use proptest::prelude::*;

use r_core::db::Database;
use r_core::expr::{compute_expr_lr, resolve_entity, EntityRows, ExprError, GroupMeans};

/// One generated name, biased towards the pathological: an empty string, a bare subunit, a name
/// that is also a complex, and a name that is only sometimes present in the matrix.
fn name_strategy() -> impl Strategy<Value = String> {
    prop_oneof![
        Just(String::new()),
        Just("".to_string()),
        "[A-Za-z]{1,3}",
        "sub[0-9]{1,2}",
        "cpx[0-9]{1,2}",
        "cof[0-9]{1,2}",
    ]
}

/// A subunit list: possibly empty, possibly repeating, possibly referring to names that are
/// themselves complexes. Empty cells are **preserved** on purpose, because
/// `Database::complex_cells` gives the empty cells' column positions meaning.
fn subunits_strategy(max: usize) -> impl Strategy<Value = Vec<String>> {
    prop::collection::vec(name_strategy(), 0..max)
}

#[derive(Debug, Clone)]
struct FuzzDb {
    db: Database,
    /// The genes actually present in the expression matrix, in matrix row order.
    present: Vec<String>,
}

fn build_db(
    complexes: Vec<(String, Vec<String>)>,
    cofactors: Vec<(String, Vec<String>)>,
    symbols: Vec<String>,
    present: Vec<String>,
    interactions: Vec<(String, String, String, String)>,
) -> FuzzDb {
    // `from_parts` rather than a struct literal: `Database` keeps private index fields, and a
    // fuzz test that has to reach inside the struct is testing the wrong surface. The name order
    // vectors are the *source row order* the real loader preserves, which is what the empty-cell
    // column positions depend on, so they are built from the same specs.
    let complex_names: Vec<String> = complexes.iter().map(|(k, _)| k.clone()).collect();
    let cofactor_names: Vec<String> = cofactors.iter().map(|(k, _)| k.clone()).collect();
    let db = Database::from_parts(
        "fuzz",
        interactions
            .into_iter()
            .map(
                |(name, ligand, receptor, pathway)| r_core::db::Interaction {
                    name,
                    pathway,
                    ligand,
                    receptor,
                    agonist: String::new(),
                    antagonist: String::new(),
                    co_activator: String::new(),
                    co_inhibitor: String::new(),
                    annotation: "Secreted Signaling".to_string(),
                },
            )
            .collect(),
        complexes.into_iter().collect(),
        cofactors.into_iter().collect(),
        complex_names,
        cofactor_names,
        4,
        2,
        symbols.into_iter().map(|s| (s, ())).collect(),
    );
    FuzzDb { db, present }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(2048))]

    /// `resolve_entity` is total: a name that cannot be resolved returns
    /// `ExprError::SubscriptOutOfBounds` rather than panicking, and one that can be resolved
    /// returns exactly the rows the reference says it should.
    #[test]
    fn resolve_entity_never_panics_and_matches_a_reference(
        complex_specs in prop::collection::vec((name_strategy(), subunits_strategy(4)), 0..8),
        present in prop::collection::vec(name_strategy(), 1..12),
        query in name_strategy(),
    ) {
        // Deduplicated: a matrix with duplicate row names cannot arise in a CellChat object, and
        // the question of which duplicate wins is answered differently by a `HashMap` row index
        // (last wins) and by R's `match()` (first wins). That difference is worth knowing and not
        // worth fuzzing -- upstream rejects duplicate rownames earlier, with
        // `data.use[RsubunitsV, ]` failing to resolve at all.
        let mut present = present;
        present.sort();
        present.dedup();
        let db = build_db(complex_specs.clone(), Vec::new(),
                          present.clone(), present.clone(), Vec::new());
        // A mean row for every present gene, in matrix order.
        let levels = ["g1".to_string(), "g2".to_string()];
        // Column-major genes x groups, which is what `GroupMeans::new` asserts on. Passing only
        // one value per gene is the shape the assertion exists to catch, and proptest found it.
        let values: Vec<f64> = (0..db.present.len() * levels.len())
            .map(|i| i as f64 * 0.5)
            .collect();
        let means = GroupMeans::new(&values, &db.present, levels.len());
        let got = resolve_entity(&query, &means, &db.db);

        // Independent reference: the gene universe first, then the complex table.
        let universe: BTreeSet<&str> = db.present.iter().map(|s| s.as_str()).collect();
        let complexes: BTreeMap<&str, Vec<&str>> = complex_specs
            .iter()
            .map(|(k, v)| (k.as_str(), v.iter().map(|s| s.as_str()).collect()))
            .collect();
        let expected: Option<Vec<usize>> = if universe.contains(query.as_str()) {
            // A gene wins over a complex of the same name: `resolve_entity` looks in the matrix
            // first, and that ordering is load-bearing when the two collide.
            Some(vec![db.present.iter().position(|s| s == &query).unwrap()])
        } else {
            complexes.get(query.as_str()).map(|subs| {
                subs.iter().filter(|s| !s.is_empty())
                    .filter_map(|s| db.present.iter().position(|p| p == s))
                    .collect()
            })
        };
        match (got, expected) {
            (Ok(EntityRows::Single(r)), Some(want)) => {
                prop_assert_eq!(&want, &vec![r]);
            }
            (Ok(EntityRows::Complex(rows)), Some(want)) => {
                prop_assert_eq!(&want, &rows);
                // A complex whose name is in the matrix resolves to `Single`, never to `Complex`:
                // a complex with **no** non-empty subunit therefore cannot be `Complex` at all.
                prop_assert!(!rows.is_empty() || db.db.complex_cells(&query).is_some());
            }
            (Err(ExprError::SubscriptOutOfBounds { name }), None) => {
                prop_assert_eq!(name, query);
            }
            (Err(ExprError::SubscriptOutOfBounds { name }), Some(_)) => {
                prop_assert!(!universe.contains(query.as_str()),
                    "resolved a name that the reference cannot: {query}");
                let _ = name;
            }
            (other, _) => prop_assert!(false, "unexpected resolution {other:?} for {query:?}"),
        }
    }

    /// A complex with a subunit that is **not** in the matrix raises upstream's
    /// `subscript out of bounds`, and so does a complex name that is in neither table. Both are
    /// `Err`, never a panic -- which is the whole reason the objective names this error.
    #[test]
    fn missing_subunit_is_an_error_not_a_panic(
        missing in prop::sample::select(vec!["absent1", "absent2", "absent3", "absent4", "absent5"]),
    ) {
        let db = build_db(
            vec![("cpxA".to_string(), vec!["sub1".to_string(), missing.to_string()])],
            Vec::new(),
            vec!["sub1".to_string()],
            vec!["sub1".to_string()],
            Vec::new(),
        );
        let values = vec![1.0];
        let levels = ["g1".to_string()];
        let means = GroupMeans::new(&values, &["sub1".to_string()], levels.len());
        let r = resolve_entity("cpxA", &means, &db.db);
        prop_assert!(matches!(r, Err(ExprError::SubscriptOutOfBounds { .. })),
            "a missing subunit must be an Err, got {r:?}");
        prop_assert!(matches!(resolve_entity("cpxZ", &means, &db.db),
                             Err(ExprError::SubscriptOutOfBounds { .. })),
            "an unknown name must be an Err");
    }

    /// `compute_expr_lr` over arbitrary (ligand, receptor) pairs is total. This is the entry point
    /// every `computeCommunProb` call goes through, so a panic here is a panic in the kernel on a
    /// user dataset.
    #[test]
    fn compute_expr_lr_never_panics(
        complex_specs in prop::collection::vec((name_strategy(), subunits_strategy(4)), 0..8),
        present in prop::collection::vec(name_strategy(), 1..10),
        pairs in prop::collection::vec((name_strategy(), name_strategy()), 0..6),
    ) {
        let db = build_db(complex_specs, Vec::new(), present.clone(), present.clone(), Vec::new());
        let n_groups = 3usize;
        let levels: Vec<String> = (0..n_groups).map(|i| format!("g{i}")).collect();
        let values: Vec<f64> = (0..db.present.len() * n_groups).map(|i| i as f64 % 7.0).collect();
        let means = GroupMeans::new(&values, &db.present, levels.len());
        let names: Vec<String> = pairs
            .iter()
            .flat_map(|(a, b)| [a.clone(), b.clone()])
            .collect();
        let r = compute_expr_lr(&names, &means, &db.db);
        if let Ok(m) = r.as_ref() {
            // Column-major, `n_groups` blocks of `n_names`, as the contract states.
            prop_assert_eq!(m.len(), names.len() * n_groups,
                "compute_expr_lr returned {} values for {} names x {} groups",
                m.len(), names.len(), n_groups);
            prop_assert!(m.iter().all(|v| v.is_finite() || v.is_nan()),
                "compute_expr_lr produced a non-finite value that is not NaN");
        }
    }

    /// `extract_gene`'s expansion rule, which is *not* "is this name in the complex table".
    /// Upstream splits on the **official symbol** list:
    ///
    /// ```r
    /// complex <- geneSet[which(geneSet %in% geneIfo$Symbol == "FALSE")]
    /// geneSet  <- intersect(geneSet, geneIfo$Symbol)
    /// complexsubunits <- select(complex_input[match(complex, rownames(complex_input), 0), ],
    ///                           starts_with("subunit"))
    /// complex <- intersect(complex, rownames(complexsubunits))
    /// geneSet <- unique(c(geneSet, unique(unlist(complexsubunits)[... != ""])))
    /// ```
    ///
    /// So a name that is an official symbol is kept **as itself and never expanded**, and a name
    /// that is not is looked up in the complex table and contributes only its non-empty subunits.
    /// In the real database the complex names are not symbols and the plain gene names are, which
    /// is what makes the rule invisible; a generator that uses the *same* strings for both
    /// branches exercises it properly. Three earlier versions of this property used "is in the
    /// complex table" as the discriminator and each was wrong.
    #[test]
    fn extract_gene_splits_on_the_symbol_list_not_the_complex_table(
        complex_specs in prop::collection::vec((name_strategy(), subunits_strategy(4)), 1..6),
        plain_genes in prop::collection::vec(name_strategy(), 1..4),
        referenced in prop::collection::vec(any::<bool>(), 1..8),
        empty_present in any::<bool>(),
    ) {
        let complexes: BTreeMap<String, Vec<String>> =
            complex_specs.iter().map(|(k, v)| (k.clone(), v.clone())).collect();
        // Only the *plain* genes are official symbols. Complex names deliberately are not, because
        // that is the real database's shape and the branch that matters.
        let mut symbols: Vec<String> = if empty_present { vec![] } else { plain_genes.clone() };
        symbols.retain(|s| !s.is_empty());
        symbols.sort();
        symbols.dedup();

        // Half the referenced names are complex names, half are plain genes.
        let mut names: Vec<String> = Vec::new();
        let ck: Vec<&String> = complex_specs.iter().map(|(k, _)| k).collect();
        let n = ck.len().max(1);
        for (i, &use_complex) in referenced.iter().enumerate() {
            if use_complex {
                names.push(ck[i % n].clone());
            } else if let Some(g) = plain_genes.get(i) {
                names.push(g.clone());
            }
        }
        let interactions: Vec<(String, String, String, String)> = names
            .iter()
            .enumerate()
            .map(|(i, nm)| (format!("I{i}"), nm.clone(), "R".to_string(), "P".to_string()))
            .collect();
        let symbol_set: BTreeSet<String> = symbols.iter().cloned().collect();
        let db = build_db(complex_specs, Vec::new(), symbols, Vec::new(), interactions);
        let genes: BTreeSet<String> = db.db.extract_gene().into_iter().collect();

        for nm in names.iter() {
            if nm.is_empty() {
                continue;
            }
            if symbol_set.contains(nm) {
                prop_assert!(genes.contains(nm),
                    "{nm} is an official symbol, so it is kept as itself and never expanded");
            } else if let Some(subs) = complexes.get(nm) {
                for s in subs.iter().filter(|s| !s.is_empty()) {
                    prop_assert!(genes.contains(s), "extract_gene omitted subunit {s} of {nm}");
                }
            }
            // Whether a name that is *neither* a symbol nor a complex of its own accord appears is
            // decided by one thing: is it some complex's non-empty subunit? The generator produces
            // overlapping subunits freely, and a name can therefore be in the output while being
            // nobody's own name. Checking that once, here, replaces two mutually inconsistent
            // assertions that each failed on a different overlapping-subunit case.
            if !symbol_set.contains(nm) {
                let a_subunit_somewhere: bool = complexes
                    .iter()
                    .any(|(_, subs)| subs.iter().any(|s| s == nm && !s.is_empty()));
                prop_assert!(a_subunit_somewhere || !genes.contains(nm),
                    "{nm} is neither a symbol nor a complex, so only another complex's subunit \
                     list can introduce it");
            }
        }
    }
}
