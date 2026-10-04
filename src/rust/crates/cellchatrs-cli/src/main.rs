//! `cellchatrs` — the CellChat inference kernel without R.
//!
//! The R package is the drop-in accelerator, but `r-core` is a plain Rust crate with no R
//! dependency, and this is the driver that makes it usable on its own: a single self-describing
//! input file in, the `Prob`/`Pval` networks and the aggregated table out, in a format that is
//! byte-comparable against what R writes for the same data.
//!
//! That last property is the point. A CLI that printed rounded numbers would be useful and
//! unverifiable; `tests/parity/check_cli.R` runs this binary against pinned upstream's
//! `computeCommunProb` on the same input and requires `identical()` on the full arrays, so "the
//! numerics are usable without R" is a checked claim rather than an aspiration.
//!
//! ```text
//! cellchatrs run       --input IN.tsv --db DIR [--out OUT.tsv] [--thresh T] [--threads N]
//! cellchatrs describe  --input IN.tsv [--db DIR]
//! cellchatrs mean      --kind KIND [--trim T] < values.txt
//! cellchatrs version
//! ```
//!
//! Values go out as `%.17g`, which round-trips every `f64` uniquely, so R's `as.numeric` on the
//! output recovers the exact bits. Inputs come in as `%a` hex floats, for the reasons in
//! [`input`].

mod input;
mod pipeline;

use std::io::Read;
use std::path::PathBuf;
use std::process::ExitCode;

const USAGE: &str = "\
cellchatrs -- the CellChat inference kernel, without R

USAGE:
  cellchatrs run       --input IN.tsv --db DIR [--out OUT.tsv] [--thresh T] [--threads N] [--hex]
  cellchatrs describe  --input IN.tsv [--db DIR]
  cellchatrs mean      --kind KIND [--trim T] VALUE...
  cellchatrs version

COMMANDS:
  run        computeProb and aggregateNet, and write the result as TSV
  describe   print the resolved shape of an input without computing anything
  mean       print one mean of the given numbers, or of stdin if none are given
  version    print the version

OPTIONS:
  --input    the input file; see tests/parity/gen_cli_fixture.R for the format
  --db       a CellChatDB export directory, as written by tests/parity/export_db.R
  --out      where to write the result; stdout if omitted
  --thresh   aggregateNet's threshold, default 0.05, matching upstream
  --threads  worker threads; defaults to the rayon pool
  --kind     triMean | truncatedMean | thresholdedMean | geometricMean
  --trim     trim fraction for the trimmed and thresholded means, default 0.1
  --hex      write the numbers as hex floats instead of %.17g; both round-trip exactly
";

/// `%.17g` is the shortest precision that round-trips every `f64` uniquely, so R's `as.numeric` on
/// this recovers the exact bit pattern. The parity gate depends on it.
fn fmt17(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    format!("{v:.17e}")
}

struct Args {
    command: String,
    input: Option<PathBuf>,
    db: Option<PathBuf>,
    out: Option<PathBuf>,
    thresh: f64,
    threads: Option<usize>,
    kind: Option<String>,
    trim: f64,
    hex: bool,
    /// Positional numbers for `cellchatrs mean`.
    values: Vec<String>,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        command: String::new(),
        input: None,
        db: None,
        out: None,
        thresh: 0.05,
        threads: None,
        kind: None,
        trim: 0.1,
        hex: false,
        values: Vec::new(),
    };
    let mut it = std::env::args().skip(1);
    a.command = it.next().unwrap_or_default();
    if a.command.is_empty() || a.command == "--help" || a.command == "-h" {
        return Err(String::new());
    }
    let rest: Vec<String> = it.collect();
    let mut i = 0;
    while i < rest.len() {
        let need = |i: usize, name: &str| -> Result<String, String> {
            rest.get(i + 1)
                .cloned()
                .ok_or_else(|| format!("{name} needs a value"))
        };
        match rest[i].as_str() {
            "--input" => {
                a.input = Some(PathBuf::from(need(i, "--input")?));
                i += 2;
            }
            "--db" => {
                a.db = Some(PathBuf::from(need(i, "--db")?));
                i += 2;
            }
            "--out" => {
                a.out = Some(PathBuf::from(need(i, "--out")?));
                i += 2;
            }
            "--thresh" => {
                a.thresh = need(i, "--thresh")?
                    .parse()
                    .map_err(|_| "--thresh needs a number".to_string())?;
                i += 2;
            }
            "--threads" => {
                a.threads = Some(
                    need(i, "--threads")?
                        .parse()
                        .map_err(|_| "--threads needs an integer".to_string())?,
                );
                i += 2;
            }
            "--kind" => {
                a.kind = Some(need(i, "--kind")?);
                i += 2;
            }
            "--hex" => {
                a.hex = true;
                i += 1;
            }
            "--trim" => {
                a.trim = need(i, "--trim")?
                    .parse()
                    .map_err(|_| "--trim needs a number".to_string())?;
                i += 2;
            }
            other if other.starts_with('-') && other.parse::<f64>().is_err() => {
                return Err(format!("unknown option {other:?}"))
            }
            other => {
                a.values.push(other.to_string());
                i += 1;
            }
        }
    }
    Ok(a)
}

/// `--hex` switches the output from `%.17g` to `%a`. Both round-trip exactly; hex is lossless by
/// construction (every value is its own bit pattern) and `%.17g` is lossless by the 17-digit
/// theorem, so the choice is readability against diffability rather than accuracy.
fn write_result(run: &pipeline::Run, inp: &input::Input, thresh: f64, hex: bool) -> String {
    let f = if hex { input::format_hex } else { fmt17 };
    let mut s = String::new();
    s.push_str("version\t1\n");
    s.push_str(&format!("groups\t{}\n", run.n_groups));
    s.push_str(&format!("lr\t{}\n", run.n_lr));
    // The input's own shape, echoed. The network's `dim` is `K x K x N` and says nothing about how
    // many genes went in, so a consumer of this file -- `tests/parity/check_cli.R` needs it to know
    // how long the `ave_expr` block is -- has no way to work it out otherwise.
    s.push_str(&format!("genes\t{}\n", inp.genes.len()));
    s.push_str(&format!("cells\t{}\n", inp.cells.len()));
    s.push_str(&format!("type_mean\t{}\n", run.type_mean_resolved));
    s.push_str(&format!(
        "dim\t{},{},{}\n",
        run.n_groups, run.n_groups, run.n_lr
    ));
    s.push_str("prob\n");
    for v in &run.prob {
        s.push_str(&f(*v));
        s.push('\n');
    }
    s.push_str("pval\n");
    for v in &run.pval {
        s.push_str(&f(*v));
        s.push('\n');
    }
    // Two averages, and they are not the same thing, so they are reported under two names.
    //
    // `ave_expr` is [`pipeline::compute_ave_expr`], which is upstream's `computeAveExpr`: it always
    // reads `data.signaling` and has no `raw.use` parameter.
    //
    // `kernel_ave` is the observed per-group aggregate the kernel actually used, which for
    // `raw.use = FALSE` comes from `data.smooth`. An earlier version reported the kernel's aggregate
    // under the name `ave_expr`, and the labels agreed for thirteen of the fourteen configurations --
    // only `raw.use = FALSE` distinguished them, and there `computeAveExpr means MISMATCH 96 values`.
    // A block whose name is wrong for one setting of one argument is a block that cannot be trusted.
    match pipeline::compute_ave_expr(inp) {
        Ok(ave) => {
            s.push_str("ave_expr\n");
            for v in &ave {
                s.push_str(&f(*v));
                s.push('\n');
            }
        }
        Err(e) => {
            s.push_str("ave_expr\tunavailable\t");
            s.push_str(&e.replace('\n', " "));
            s.push('\n');
        }
    }
    s.push_str("kernel_ave\n");
    for v in &run.avg {
        s.push_str(&f(*v));
        s.push('\n');
    }
    let (count, weight) =
        pipeline::aggregate_net(&run.prob, &run.pval, run.n_groups, run.n_lr, thresh);
    s.push_str("aggregate_count\n");
    for row in &count {
        s.push_str(&row.iter().map(|v| f(*v)).collect::<Vec<_>>().join("\t"));
        s.push('\n');
    }
    s.push_str("aggregate_weight\n");
    for row in &weight {
        s.push_str(&row.iter().map(|v| f(*v)).collect::<Vec<_>>().join("\t"));
        s.push('\n');
    }
    s.push_str("dimnames_source\n");
    s.push_str(&inp.groups.join("\t"));
    s.push('\n');
    s.push_str("dimnames_target\n");
    s.push_str(&inp.groups.join("\t"));
    s.push('\n');
    s.push_str("dimnames_interaction\n");
    s.push_str(
        &inp.lr
            .iter()
            .map(|p| p.label.clone())
            .collect::<Vec<_>>()
            .join("\t"),
    );
    s.push('\n');
    s
}

fn run() -> Result<(), String> {
    let a = parse_args().map_err(|e| {
        if e.is_empty() {
            USAGE.to_string()
        } else {
            format!("{e}\n\n{USAGE}")
        }
    })?;

    match a.command.as_str() {
        "version" | "--version" => {
            println!("cellchatrs {}", env!("CARGO_PKG_VERSION"));
            println!("r-core {}", r_core::VERSION);
            return Ok(());
        }
        "mean" => {
            let kind = a
                .kind
                .clone()
                .ok_or("--kind needs a value; see --help".to_string())?;
            // Values may be given as arguments or on stdin. The arguments come first because
            // `system2` does not forward stdin: a test harness (and `check_cli.R`) runs this with
            // `R --vanilla -f script.R`, where the child inherits an already-consumed script as its
            // stdin and silently reads zero numbers, which surfaces as a `NaN` rather than an error.
            // Reading arguments makes the subcommand usable from a shell pipeline *and* testable.
            let text = if a.values.is_empty() {
                let mut buf = String::new();
                std::io::stdin()
                    .read_to_string(&mut buf)
                    .map_err(|e| format!("cannot read stdin: {e}"))?;
                buf
            } else {
                a.values.join(" ")
            };
            let x: Vec<f64> = text
                .split_whitespace()
                .map(|t| {
                    t.parse::<f64>()
                        .map_err(|_| format!("{t:?} is not a number"))
                })
                .collect::<Result<_, _>>()?;
            if x.is_empty() {
                return Err("mean needs at least one value, as arguments or on stdin".into());
            }
            let out = pipeline::means(&x, &kind, a.trim)?;
            for v in out {
                println!("{}", fmt17(v));
            }
            return Ok(());
        }
        _ => {}
    }

    let path = a
        .input
        .clone()
        .ok_or("run and describe need --input; see --help".to_string())?;
    let inp = input::read(&path).map_err(|e| format!("{path:?}: {e}"))?;

    if let Some(t) = a.threads {
        // `build_global` is ignored if the pool already exists, which is the right behaviour: a
        // second call in one process is a no-op rather than a panic, and the CLI is single-shot.
        rayon::ThreadPoolBuilder::new()
            .num_threads(t)
            .build_global()
            .map_err(|e| format!("cannot set the thread count to {t}: {e}"))?;
    }

    if a.command == "describe" {
        println!("groups\t{}", inp.groups.len());
        println!("genes\t{}", inp.genes.len());
        println!("cells\t{}", inp.cells.len());
        println!("lr\t{}", inp.lr.len());
        println!("type\t{}", inp.cfg.type_mean);
        println!("trim\t{}", fmt17(inp.cfg.trim));
        println!("population_size\t{}", inp.cfg.population_size);
        println!("raw_use\t{}", inp.cfg.raw_use);
        println!("nboot\t{}", inp.cfg.nboot);
        println!("seed\t{}", inp.cfg.seed);
        println!("Kh\t{}", fmt17(inp.cfg.kh));
        println!("n\t{}", fmt17(inp.cfg.n));
        for (i, l) in inp.groups.iter().enumerate() {
            let n = inp.cell_group.iter().filter(|&&g| g == i).count();
            println!("group[{i}]\t{l}\t{n} cells");
        }
        // The observed per-group means, so `describe` can show what the kernel is about to read
        // without a full run. `computeAveExpr` is part of the required numeric surface, and a
        // describe that stops at the shape leaves it unreachable from the CLI.
        match pipeline::compute_ave_expr(&inp) {
            Ok(avg) => {
                println!("ave_expr\t{} values", avg.len());
                let n_genes = inp.genes.len();
                let n_groups = inp.groups.len();
                for g in 0..n_genes {
                    let row: Vec<String> =
                        (0..n_groups).map(|c| fmt17(avg[c * n_genes + g])).collect();
                    println!("ave_expr[{}]\t{}", inp.genes[g], row.join("\t"));
                }
            }
            Err(e) => println!("ave_expr\tunavailable: {e}"),
        }
        if let Some(db) = &a.db {
            let d = r_core::db::Database::load(db).map_err(|e| format!("{db:?}: {e}"))?;
            println!("db_species\t{}", d.species);
            println!("db_interactions\t{}", d.interactions.len());
            println!("db_complexes\t{}", d.complexes.len());
            println!("db_cofactors\t{}", d.cofactors.len());
        }
        return Ok(());
    }

    if a.command != "run" {
        return Err(format!("unknown command {:?}\n\n{USAGE}", a.command));
    }

    let db_path = a
        .db
        .clone()
        .ok_or("run needs --db; the interaction table alone is not enough to resolve complexes")?;
    let db = r_core::db::Database::load(&db_path).map_err(|e| format!("{db_path:?}: {e}"))?;
    let run = pipeline::compute_commun_prob(&inp, &db)?;
    let text = write_result(&run, &inp, a.thresh, a.hex);
    match &a.out {
        Some(p) => std::fs::write(p, text).map_err(|e| format!("cannot write {p:?}: {e}")),
        None => {
            print!("{text}");
            Ok(())
        }
    }
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(msg) => {
            eprintln!("cellchatrs: {msg}");
            ExitCode::FAILURE
        }
    }
}
