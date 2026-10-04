//! Parity of [`r_core::db`] against the pinned `CellChatDB`.
//!
//! Fixtures under `tests/fixtures/db_{human,mouse}/` are produced by
//! `tests/parity/export_db.R`, which also evaluates **upstream's own R code**
//! (`extractGene` / `extractGeneSubset`, verbatim from `R/database.R`) and writes the
//! result to `extract_gene*.txt`. The Rust port is compared against that
//! order-for-order — not merely set-for-set, because `identifyOverExpressedGenes` and
//! `computeAveExpr(features = ...)` observe the order.
//!
//! Manifests record the source `.rda`'s MD5 and the upstream commit, so a stale or
//! doctored fixture fails loudly rather than silently.

use r_core::db::{Database, ANNOTATION_LEVELS};
use std::path::{Path, PathBuf};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/fixtures")
}

fn read_lines(p: &Path) -> Vec<String> {
    std::fs::read_to_string(p)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", p.display()))
        .lines()
        .map(|s| s.to_string())
        .collect()
}

fn manifest(dir: &Path) -> Vec<(String, String)> {
    read_lines(&dir.join("MANIFEST.tsv"))
        .iter()
        .filter_map(|l| {
            let l = l.trim();
            if l.starts_with('#') {
                None
            } else {
                l.split_once('\t')
                    .map(|(a, b)| (a.to_string(), b.to_string()))
            }
        })
        .collect()
}

const SPECIES: [&str; 2] = ["human", "mouse"];

fn load(species: &str) -> (Database, PathBuf) {
    let dir = fixtures().join(format!("db_{species}"));
    let db = Database::load(&dir).unwrap_or_else(|e| panic!("{species}: {e}"));
    db.verify_manifest(&dir)
        .unwrap_or_else(|e| panic!("{species}: {e}"));
    (db, dir)
}

#[test]
fn both_species_load_and_match_the_manifest_counts() {
    for sp in SPECIES {
        let (db, dir) = load(sp);
        let m = manifest(&dir);
        let get = |k: &str| {
            m.iter()
                .find(|(a, _)| a == k)
                .unwrap_or_else(|| panic!("{sp}: manifest lacks {k}"))
                .1
                .clone()
        };
        assert_eq!(db.n_lr().to_string(), get("n_interactions"), "{sp}: nLR");
        assert_eq!(
            db.complexes.len().to_string(),
            get("n_complexes"),
            "{sp}: complexes"
        );
        assert_eq!(
            db.cofactors.len().to_string(),
            get("n_cofactors"),
            "{sp}: cofactors"
        );
        assert_eq!(
            db.symbols.len().to_string(),
            get("n_symbols"),
            "{sp}: symbols"
        );
        assert_eq!(db.nlr1.to_string(), get("nlr1"), "{sp}: nLR1");
    }
}

#[test]
fn extract_gene_matches_upstream_r_order_for_order() {
    for sp in SPECIES {
        let (db, dir) = load(sp);
        let want = read_lines(&dir.join("extract_gene.txt"));
        let got = db.extract_gene();
        assert_eq!(
            got.len(),
            want.len(),
            "{sp}: length {} vs R's {}",
            got.len(),
            want.len()
        );
        for (i, (g, w)) in got.iter().zip(&want).enumerate() {
            assert_eq!(g, w, "{sp}: extract_gene differs at index {i}");
        }
    }
}

/// The default `subsetDB` drops `Non-protein Signaling`, which changes both the
/// interaction count and the gene set. Getting this wrong changes `data.signaling`,
/// hence every group mean, hence the whole network.
#[test]
fn extract_gene_matches_for_the_three_annotation_subset() {
    for sp in SPECIES {
        let (db, dir) = load(sp);
        let want = read_lines(&dir.join("extract_gene_three.txt"));
        let three: Vec<String> = ANNOTATION_LEVELS
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != 2)
            .map(|(_, l)| l.to_string())
            .collect();

        let mut sub = db.clone();
        sub.interactions.retain(|i| three.contains(&i.annotation));
        let got = sub.extract_gene();
        assert_eq!(got, want, "{sp}: 3-annotation extract_gene differs");
        // Sanity: the 4-annotation set is a strict superset here, and strictly bigger.
        assert!(want.len() < db.extract_gene().len());
    }
}

/// `extractGeneSubset`'s "is this a complex?" test is `x %in% symbols == "FALSE"`, which
/// works only because R coerces the logical to character. Reading it literally drops
/// every complex and its subunits; assert they survive.
#[test]
fn complex_subunits_are_reachable_through_extract_gene_subset() {
    let (db, _) = load("human");
    // TGFbR1_R2 is a receptor complex in the human DB.
    let out = db.extract_gene_subset(&["TGFbR1_R2".to_string(), "TGFB1".to_string()]);
    assert!(
        out.contains(&"TGFB1".to_string()),
        "the plain gene must be kept"
    );
    let subunits = db
        .complex_subunits("TGFbR1_R2")
        .expect("TGFbR1_R2 must be in the complex table");
    assert!(!subunits.is_empty(), "complex must have subunits");
    for s in subunits {
        assert!(
            out.iter().any(|g| g == s),
            "subunit {s} of TGFbR1_R2 was dropped by extract_gene_subset"
        );
    }
}

#[test]
fn a_non_official_name_with_no_complex_row_yields_nothing() {
    let (db, _) = load("human");
    let out = db.extract_gene_subset(&["NOT_A_REAL_COMPLEX_OR_GENE".to_string()]);
    assert!(
        out.is_empty(),
        "an unknown non-symbol must contribute nothing, got {out:?}"
    );
}

#[test]
fn annotation_order_is_a_contiguous_stable_partition() {
    for sp in SPECIES {
        let (db, _) = load(sp);
        // Exactly the ordering `computeCommunProb` depends on for `nLR1`.
        let order: Vec<usize> = ANNOTATION_LEVELS
            .iter()
            .enumerate()
            .flat_map(|(li, l)| {
                db.interactions
                    .iter()
                    .enumerate()
                    .filter(move |(_, i)| i.annotation == *l)
                    .map(move |(i, _)| (li, i))
            })
            .map(|(_, i)| i)
            .collect();
        let actual: Vec<usize> = (0..db.n_lr()).collect();
        assert_eq!(
            order, actual,
            "{sp}: interactions are not in annotation order"
        );

        // And each level occupies a contiguous run, so the diffusion/contact split is a
        // single boundary.
        let mut runs: Vec<&str> = Vec::new();
        for i in &db.interactions {
            if runs.last() != Some(&i.annotation.as_str()) {
                runs.push(&i.annotation);
            }
        }
        assert_eq!(
            runs.len(),
            ANNOTATION_LEVELS.len(),
            "{sp}: expected {} contiguous annotation runs, got {runs:?}",
            ANNOTATION_LEVELS.len()
        );
        // nLR1 is the last index of the first run.
        assert_eq!(db.nlr1, db.diffusion_mediated().count(), "{sp}: nLR1");
    }
}

#[test]
fn every_interaction_ligand_or_receptor_is_a_symbol_or_a_complex() {
    // Otherwise `computeExpr_LR` has no path for it: a name that is neither an official
    // symbol nor a complex row would be routed to `computeExpr_complex` and then
    // subscript out of bounds (R-ism 11).
    for sp in SPECIES {
        let (db, _) = load(sp);
        for it in db.interactions.iter() {
            for (role, name) in [("ligand", &it.ligand), ("receptor", &it.receptor)] {
                assert!(
                    db.is_symbol(name) || db.complex_subunits(name).is_some(),
                    "{sp}: interaction {} {role} {name:?} is neither a symbol nor a complex",
                    it.name
                );
                if let Some(s) = db.complex_subunits(name) {
                    for sub in s {
                        assert!(
                            db.is_symbol(sub),
                            "{sp}: interaction {} {role} complex {name:?} has non-symbol subunit {sub:?}",
                            it.name
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn cofactor_names_all_resolve() {
    for sp in SPECIES {
        let (db, _) = load(sp);
        for it in &db.interactions {
            for (role, name) in [
                ("agonist", &it.agonist),
                ("antagonist", &it.antagonist),
                ("co_activator", &it.co_activator),
                ("co_inhibitor", &it.co_inhibitor),
            ] {
                if name.is_empty() {
                    continue;
                }
                assert!(
                    db.cofactor_subunits(name).is_some(),
                    "{sp}: interaction {} {role} {name:?} is not in the cofactor table",
                    it.name
                );
            }
        }
    }
}

#[test]
fn complex_and_cofactor_subunits_are_themselves_symbols() {
    for sp in SPECIES {
        let (db, _) = load(sp);
        for (name, cells) in db.complexes.iter().chain(db.cofactors.iter()) {
            for s in cells.iter().filter(|s| !s.is_empty()) {
                assert!(
                    db.is_symbol(s),
                    "{sp}: subunit {s:?} of {name:?} is not an official symbol"
                );
            }
        }
    }
}

#[test]
fn a_tampered_export_is_rejected() {
    // The staleness check has to actually bite, or it is decoration.
    let (db, dir) = load("human");
    let tampered = fixtures().join("db_human_tampered");
    let _ = std::fs::remove_dir_all(&tampered);
    for f in [
        "interactions.tsv",
        "complexes.tsv",
        "cofactors.tsv",
        "symbols.txt",
        "MANIFEST.tsv",
    ] {
        let text = std::fs::read_to_string(dir.join(f)).unwrap();
        std::fs::create_dir_all(&tampered).unwrap();
        std::fs::write(tampered.join(f), text).unwrap();
    }
    // Flip one gene symbol; the manifest's MD5 must notice.
    let p = tampered.join("symbols.txt");
    let text = std::fs::read_to_string(&p).unwrap();
    std::fs::write(&p, text.replace("TGFB1", "TGFZ9")).unwrap();
    let err = Database::load(&tampered).unwrap_err();
    assert!(
        format!("{err}").contains("stale"),
        "expected a staleness error, got: {err}"
    );
    let _ = std::fs::remove_dir_all(&tampered);
    let _ = db;
}
