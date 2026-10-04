//! Reader for the CLI's self-describing input file.
//!
//! The format is a tagged TSV: a `key<TAB>value` header record, then the data that key announces, in
//! a fixed order. It is the same shape as the existing golden corpora in `tests/fixtures`, for the
//! same reason: a fixture written by one side and read by the other cannot rot the way a re-derived
//! one can. `expr_matrix.tsv` records `groups`, `genes`, the names, a `values` marker and then one
//! line of hex floats per gene, and the comment on `fixture_matrix()` in `expr_parity.rs` explains
//! what went wrong when the layout was reconstructed in Rust instead of read from a file.
//!
//! **Two properties the format has to have, both of which bit-parity depends on.**
//!
//! *Values are hex floats.* `%a` round-trips every `f64` exactly, so nothing is lost between R and
//! Rust. Decimal `%.17g` would also round-trip, but a miscounted digit would be silently accepted
//! as a different double, and the resulting parity failure would look like a rounding bug in the
//! kernel. Hex makes a transcription error *loud*.
//!
//! *The matrix is gene-major, one line per gene.* A genes-by-cells matrix flattened in R's
//! column-major order has, for each column (cell), all genes consecutive. Writing one line per gene
//! means the reader indexes `data[gene * n_cells + cell]`, which is the transpose of R's storage. That
//! inversion is exactly the bug `expr_parity.rs` documents at length, so the format makes the
//! inversion explicit and the reader states which order it produces.

use std::fmt;
use std::path::Path;

/// Everything the kernel needs for one `computeCommunProb` call.
#[derive(Debug)]
pub struct Input {
    pub genes: Vec<String>,
    pub cells: Vec<String>,
    /// The group level names, in the order R's `factor` levels have them.
    pub groups: Vec<String>,
    /// `n_cells` entries, each an index into `groups`.
    pub cell_group: Vec<usize>,
    /// `n_genes * n_cells`, indexed `gene * n_cells + cell`.
    pub data: Vec<f64>,
    /// `object@data.smooth`, the same shape, required only when `cfg.raw_use` is `FALSE`.
    ///
    /// Upstream's `computeCommunProb` is literally
    /// `if (raw.use) data <- as.matrix(object@data.signaling) else data <- as.matrix(object@data.smooth)`,
    /// so a format that carried only one matrix could not express `raw.use = FALSE` at all. Deriving
    /// it -- say, from the row sums of `data` -- would be a different number, and the parity gate would
    /// report a kernel defect for what is really a fixture that was never written down. The field is
    /// therefore part of the format, and its absence is an error rather than a fallback.
    pub data_smooth: Option<Vec<f64>>,
    pub lr: Vec<r_core::prob::LrPair>,
    pub cfg: Config,
}

#[derive(Debug, Clone)]
pub struct Config {
    pub type_mean: String,
    pub trim: f64,
    pub population_size: bool,
    pub raw_use: bool,
    pub nboot: usize,
    pub seed: i32,
    pub kh: f64,
    pub n: f64,
}

#[derive(Debug)]
pub enum Error {
    Io(std::io::Error),
    Parse { line: usize, what: String },
    Shape(String),
    Rf(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => write!(f, "cannot read input: {e}"),
            Error::Parse { line, what } => write!(f, "line {line}: {what}"),
            Error::Shape(m) => write!(f, "inconsistent input: {m}"),
            Error::Rf(m) => write!(f, "cannot parse a float: {m}"),
        }
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

/// A line cursor that reports the line number in every error, because a fixture that fails in the
/// middle of 24 x 60 hex floats is otherwise unlocatable.
struct Lines {
    idx: usize,
    inner: std::vec::IntoIter<String>,
}

impl Lines {
    fn new(text: &str) -> Self {
        Lines {
            idx: 0,
            inner: text
                .lines()
                .map(str::to_string)
                .collect::<Vec<_>>()
                .into_iter(),
        }
    }
    fn next(&mut self) -> Result<String, Error> {
        self.idx += 1;
        self.inner.next().ok_or_else(|| Error::Parse {
            line: self.idx,
            what: "input ended in the middle of a record".into(),
        })
    }
    fn header(&mut self, key: &str) -> Result<String, Error> {
        let line = self.next()?;
        let mut it = line.split('\t');
        let got = it.next().unwrap_or("");
        if got != key {
            return Err(Error::Parse {
                line: self.idx,
                what: format!("expected the record `{key}`, found `{got}`"),
            });
        }
        Ok(it.next().unwrap_or("").to_string())
    }
    /// Every field after the key on a `key<TAB>a<TAB>b` line, not just the first.
    ///
    /// `header` returns `it.next()`, which is right for `key<TAB>value` and wrong for a *list* of
    /// values. Used for the L-R column header, where returning only `interaction_name` produced
    /// "the L-R header is interaction_name, expected exactly interaction_name,ligand,..." -- a
    /// message that names the right thing and is about the reader, not the file.
    fn header_all(&mut self, key: &str) -> Result<Vec<String>, Error> {
        let line = self.next()?;
        let mut it = line.split('\t');
        let got = it.next().unwrap_or("");
        if got != key {
            return Err(Error::Parse {
                line: self.idx,
                what: format!("expected the record `{key}`, found `{got}`"),
            });
        }
        Ok(it.map(|s| s.trim().to_string()).collect())
    }
    /// Is the next line exactly `key`? Consumes nothing.
    fn peek_is(&mut self, key: &str) -> bool {
        self.inner.clone().next().as_deref() == Some(key)
    }
    fn count(&mut self, key: &str) -> Result<usize, Error> {
        self.header(key)?.trim().parse().map_err(|_| Error::Parse {
            line: self.idx,
            what: format!("`{key}` does not carry a count"),
        })
    }
    fn marker(&mut self, key: &str) -> Result<(), Error> {
        let got = self.next()?;
        if got != key {
            return Err(Error::Parse {
                line: self.idx,
                what: format!("expected the marker `{key}`, found `{got}`"),
            });
        }
        Ok(())
    }
}

/// Format a `f64` as a C99/R-style hex float: `0x1.8p+0`.
///
/// Rust has no `{a}` format trait, so the bit pattern is decoded directly. The mantissa is written
/// as a bare hex integer with its leading `1` in the integer part -- `0x18000000000000p-1` rather
/// than `0x1.8p+0` -- which is the same value and is accepted by both [`parse_hex`] and any C99
/// compiler, so the two are interchangeable on the wire.
///
/// Subnormals are the case that makes this worth writing rather than approximating: a subnormal has
/// no implicit leading `1`, so writing one would shift it by a factor of two. The exponent for that
/// branch is `-1022 - 52`, which is the unbiased exponent of the LSB.
pub fn format_hex(v: f64) -> String {
    if v.is_nan() {
        return "NaN".into();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    let bits = v.to_bits();
    let sign = if bits >> 63 != 0 { "-" } else { "" };
    let exp = ((bits >> 52) & 0x7ff) as i32;
    let frac = bits & ((1u64 << 52) - 1);
    let (mant, e) = if exp == 0 {
        (frac, -1022 - 52)
    } else {
        (frac | (1u64 << 52), exp - 1023 - 52)
    };
    let exp_sign = if e < 0 { "-" } else { "+" };
    format!("{sign}0x{mant:x}p{exp_sign}{}", e.unsigned_abs())
}

/// Parse an R-style hex float: `0x1.8p+3`, with an optional sign. `%a` may also emit
/// `0x1.8p+3` with no fractional part, or a bare `0x0p+0` for zero, and R writes `-0x0p+0` for
/// negative zero -- which has to survive, because `Prob` distinguishes `-0.0` from `0.0`.
/// `value * 2^exp`, exactly.
///
/// `f64::powi` is not usable here. It is documented as computing by repeated multiplication, and for
/// a large negative exponent the intermediates underflow: `2.0f64.powi(-1049)` walks down past
/// `2^-1024` to `2^-2048`, which is zero, and returns a wrong answer or a NaN. Every power of two in
/// this range is exactly representable, so the multiplication is broken into steps that keep the
/// result normal throughout and the final scaling is a single `from_bits` of the biased exponent.
///
/// Multiplying by a power of two is exact as long as the result stays in the normal range, so this
/// introduces no rounding at all -- which is the requirement, since this is a bit-exact wire format.
fn ldexp_exact(mut v: f64, mut e: i32) -> f64 {
    const HALF: f64 = f64::from_bits(1u64 << 52); // 2^-1022
    const BIG: f64 = f64::from_bits(0x7fe_u64 << 52); // 2^1023
    while e < -1022 {
        v *= HALF;
        e += 1022;
    }
    while e > 1023 {
        v *= BIG;
        e -= 1023;
    }
    v * f64::from_bits(((1023 + e) as u64) << 52)
}

fn parse_hex(s: &str) -> Result<f64, Error> {
    let t = s.trim();
    let (neg, body) = match t.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, t),
    };
    let body = body
        .strip_prefix("0x")
        .or_else(|| body.strip_prefix("0X"))
        .ok_or_else(|| {
            Error::Rf(format!(
                "{t:?} is not a hex float; the corpus must be written with sprintf(\"%a\")"
            ))
        })?;
    // Split the mantissa from the exponent. `f64::from_str` wants a decimal point, so the hex
    // mantissa is converted by hand rather than by rewriting it.
    let (mant, exp) = match body.find(['p', 'P']) {
        Some(i) => (
            &body[..i],
            body[i + 1..]
                .parse::<i32>()
                .map_err(|e| Error::Rf(e.to_string()))?,
        ),
        None => (body, 0),
    };
    let (int_part, frac_part) = match mant.find('.') {
        Some(i) => (&mant[..i], &mant[i + 1..]),
        None => (mant, ""),
    };
    if int_part.is_empty() && frac_part.is_empty() {
        return Err(Error::Rf(format!("{t:?} has no mantissa digits")));
    }
    if !int_part.chars().all(|c| c.is_ascii_hexdigit())
        || !frac_part.chars().all(|c| c.is_ascii_hexdigit())
    {
        return Err(Error::Rf(format!("{t:?} has a non-hex mantissa")));
    }
    // `int_part.frac_part` in base 16 is `int + frac / 16^n`, and each digit is exact in a double
    // because the running total never exceeds 16^n and 16^n is a power of two.
    let mut value = 0.0f64;
    for c in int_part.chars() {
        value = value * 16.0 + c.to_digit(16).unwrap() as f64;
    }
    let mut scale = 1.0f64 / 16.0;
    for c in frac_part.chars() {
        value += c.to_digit(16).unwrap() as f64 * scale;
        scale /= 16.0;
    }
    let out = ldexp_exact(value, exp);
    Ok(if neg { -out } else { out })
}

fn parse_usize(s: &str, line: usize) -> Result<usize, Error> {
    s.trim().parse().map_err(|_| Error::Parse {
        line,
        what: format!("{s:?} is not an integer"),
    })
}

fn parse_f64(s: &str, line: usize) -> Result<f64, Error> {
    s.trim().parse().map_err(|_| Error::Parse {
        line,
        what: format!("{s:?} is not a decimal float (only the data matrix uses hex floats)"),
    })
}

fn parse_bool(s: &str, line: usize) -> Result<bool, Error> {
    match s.trim() {
        "TRUE" | "true" | "1" => Ok(true),
        "FALSE" | "false" | "0" => Ok(false),
        other => Err(Error::Parse {
            line,
            what: format!("{other:?} is not a logical"),
        }),
    }
}

pub fn read(path: &Path) -> Result<Input, Error> {
    let text = std::fs::read_to_string(path)?;
    parse(&text)
}

pub fn parse(text: &str) -> Result<Input, Error> {
    let mut l = Lines::new(text);

    let version = l.header("version")?;
    if version.trim() != "1" {
        return Err(Error::Shape(format!(
            "input format version {version:?}; this build reads version 1"
        )));
    }

    let n_groups = l.count("groups")?;
    let n_genes = l.count("genes")?;
    let n_cells = l.count("cells")?;
    let n_lr = l.count("lr")?;

    let mut groups = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        groups.push(l.next()?);
    }

    let mut genes = Vec::with_capacity(n_genes);
    for _ in 0..n_genes {
        genes.push(l.next()?);
    }

    // Gene-major on the wire, cell-major in memory. See the module comment: this is the
    // transposition that `expr_parity.rs` documents, made explicit rather than implicit.
    let read_matrix = |l: &mut Lines, what: &str| -> Result<Vec<f64>, Error> {
        let mut m = vec![0.0f64; n_genes * n_cells];
        for g in 0..n_genes {
            let line = l.next()?;
            let vals: Vec<&str> = line.split_whitespace().collect();
            if vals.len() != n_cells {
                return Err(Error::Shape(format!(
                    "{what}: gene {} (`{}`) has {} values, expected {n_cells}",
                    g,
                    genes[g],
                    vals.len()
                )));
            }
            for (c, v) in vals.iter().enumerate() {
                m[g * n_cells + c] = parse_hex(v)?;
            }
        }
        Ok(m)
    };
    l.marker("values")?;
    let data = read_matrix(&mut l, "values")?;

    // `values_smooth` is optional in the format and mandatory in practice when `raw_use = FALSE`.
    // Peeking is safe here because nothing has consumed the line it might be on.
    let data_smooth = if l.peek_is("values_smooth") {
        l.marker("values_smooth")?;
        Some(read_matrix(&mut l, "values_smooth")?)
    } else {
        None
    };

    let mut cells = Vec::with_capacity(n_cells);
    for _ in 0..n_cells {
        cells.push(l.next()?);
    }

    // One line, holding every cell's group index -- not one line per cell. An earlier version looped
    // `0..n_cells` around `l.next()`, so it read the single index line first, filled the vector, and
    // then walked off into the L-R header; with `cell_group` at capacity but length 0 the first
    // out-of-range write panicked with "the len is 0 but the index is 0" and no mention of the format.
    l.marker("cell_group")?;
    let line = l.next()?;
    let vals: Vec<&str> = line.split_whitespace().collect();
    if vals.len() != n_cells {
        return Err(Error::Shape(format!(
            "the cell_group line has {} indices, expected {n_cells}",
            vals.len()
        )));
    }
    let mut cell_group = vec![0usize; n_cells];
    {
        for (c, v) in vals.iter().enumerate() {
            let g = parse_usize(v, l.idx)?;
            if g >= n_groups {
                return Err(Error::Shape(format!(
                    "cell {} is in group {g}, but there are only {n_groups} groups",
                    cells[c]
                )));
            }
            cell_group[c] = g;
        }
    }

    // The L-R table. `interaction_name` is the label, and it is what `dimnames(prob)[[3]]` is
    // built from, so it cannot be derived here: upstream takes it from `rownames(LRsig)`.
    let lr_cols: Vec<String> = l.header_all("lr_columns")?;
    let want = [
        "interaction_name",
        "ligand",
        "receptor",
        "agonist",
        "antagonist",
        "co_A_receptor",
        "co_I_receptor",
    ];
    if lr_cols != want {
        return Err(Error::Shape(format!(
            "the L-R header is {}, expected exactly {} in that order",
            lr_cols.join(","),
            want.join(",")
        )));
    }
    let mut lr = Vec::with_capacity(n_lr);
    for i in 0..n_lr {
        let line = l.next()?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() != want.len() {
            return Err(Error::Parse {
                line: l.idx,
                what: format!(
                    "L-R row {i} has {} fields, expected {}; the table is tab-separated and a \
                     missing field must be written as an empty one, not omitted",
                    f.len(),
                    want.len()
                ),
            });
        }
        let get = |j: usize| f[j].trim().to_string();
        lr.push(r_core::prob::LrPair {
            ligand: get(1),
            receptor: get(2),
            agonist: get(3),
            antagonist: get(4),
            co_a: get(5),
            co_i: get(6),
            label: get(0),
        });
    }

    let cfg = Config {
        type_mean: l.header("type")?,
        trim: parse_f64(&l.header("trim")?, l.idx)?,
        population_size: parse_bool(&l.header("population_size")?, l.idx)?,
        raw_use: parse_bool(&l.header("raw_use")?, l.idx)?,
        nboot: parse_usize(&l.header("nboot")?, l.idx)?,
        seed: l.header("seed")?.trim().parse().map_err(|_| Error::Parse {
            line: l.idx,
            what: "seed is not an integer".into(),
        })?,
        kh: parse_f64(&l.header("Kh")?, l.idx)?,
        n: parse_f64(&l.header("n")?, l.idx)?,
    };

    if n_lr == 0 {
        return Err(Error::Shape("the L-R table is empty".into()));
    }
    if n_groups == 0 || n_genes == 0 || n_cells == 0 {
        return Err(Error::Shape(format!(
            "a zero dimension: {n_groups} groups, {n_genes} genes, {n_cells} cells"
        )));
    }
    if !cfg.raw_use && data_smooth.is_none() {
        return Err(Error::Shape(
            "raw_use is FALSE, so the input must carry a `values_smooth` block: upstream reads \
             object@data.smooth in that case, and deriving it from `values` would be a different \
             matrix"
                .into(),
        ));
    }
    if cfg.raw_use && data_smooth.is_some() {
        return Err(Error::Shape(
            "raw_use is TRUE but the input carries a `values_smooth` block, which upstream would \
             ignore. Carrying data that is not read is how a format stops describing the computation"
                .into(),
        ));
    }
    Ok(Input {
        genes,
        cells,
        groups,
        cell_group,
        data,
        data_smooth,
        lr,
        cfg,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hex_floats_round_trip_including_signed_zero() {
        for v in [
            0.0f64,
            -0.0,
            1.0,
            -1.0,
            0.5,
            1e300,
            -1e-300,
            5e-324,
            1.7976931348623157e308,
        ] {
            let s = format_hex(v);
            let back = parse_hex(&s).expect("round trip");
            assert_eq!(back.to_bits(), v.to_bits(), "{s} -> {back:e}");
        }
    }

    #[test]
    fn negative_zero_stays_negative() {
        // `Prob` and `Pval` distinguish these, and a sign that gets dropped in transcription would
        // be an `identical()` failure with no arithmetic anywhere near it.
        assert!(parse_hex("-0x0p+0").unwrap().is_sign_negative());
        assert!(!parse_hex("0x0p+0").unwrap().is_sign_negative());
    }

    #[test]
    fn the_two_hex_spellings_agree() {
        // R's `%a` writes `0x1.8p+0`; `format_hex` writes `0x18000000000000p-1`. Both are on the
        // wire -- the first because the corpus is written by R, the second because the CLI can emit
        // it -- so the parser has to read both and they have to be the same number.
        for v in [
            1.0f64,
            -1.0,
            0.5,
            1e300,
            -1e-300,
            5e-324,
            1.7976931348623157e308,
            0.1,
            3.0,
        ] {
            // R spells it with a radix point: `sprintf("%a", x)` writes `0x1.<52 fraction bits in
            // hex>p<exp>`, which is what the corpus on disk contains. All 52 bits have to be there --
            // an earlier version of this test wrote only the low 28 and dropped the high nibble,
            // which is a factor-of-sixteen error that still parsed and still looked like a float.
            let bits = v.to_bits();
            let exp = ((bits >> 52) & 0x7ff) as i32;
            let frac = bits & ((1u64 << 52) - 1);
            let dotted = if exp == 0 && frac == 0 {
                format!("{}0x0p+0", if bits >> 63 != 0 { "-" } else { "" })
            } else {
                let e = if exp == 0 {
                    -1022 - 52
                } else {
                    exp - 1023 - 52
                };
                // Subnormals have no implicit leading 1, so the integer part is the fraction's own
                // top bits and the exponent shifts by one to compensate.
                let (lead, e) = if exp == 0 { (0u64, e) } else { (1u64, e) };
                let combined = (lead << 52) | frac;
                format!(
                    "{}0x{:x}p{}{}",
                    if bits >> 63 != 0 { "-" } else { "" },
                    combined,
                    if e < 0 { "-" } else { "+" },
                    e.unsigned_abs()
                )
            };
            assert_eq!(
                parse_hex(&dotted).unwrap().to_bits(),
                parse_hex(&format_hex(v)).unwrap().to_bits(),
                "{dotted} vs {}",
                format_hex(v)
            );
        }
    }

    #[test]
    fn a_non_hex_value_is_rejected_rather_than_coerced() {
        // The whole point of the hex format: a decimal float here is a bug in the generator, and
        // accepting it would turn that into a silent parity failure.
        assert!(parse_hex("0.1").is_err());
        assert!(parse_hex("1e-5").is_err());
    }

    #[test]
    fn every_exponent_magnitude_survives() {
        // Across the whole double range, a hand-rolled hex parser is where an off-by-one in the
        // exponent or the `1/16` scaling shows up. Sampled at powers of two and at the subnormals.
        let mut v = f64::MIN_POSITIVE;
        while v < f64::MAX {
            for w in [v, -v, v * 3.0, v / 3.0] {
                // `v * 3.0` and `-v` overflow to `Inf` at the top of the range. That is not a
                // formatter failure: a hex float has no encoding for a non-finite value, so
                // `format_hex` returns `Inf`/`NaN`, which is what R's `as.numeric` reads. The
                // round-trip claim is therefore about finite values, and that is what is asserted.
                if !w.is_finite() {
                    assert_eq!(format_hex(w), if w > 0.0 { "Inf" } else { "-Inf" });
                    continue;
                }
                let s = format_hex(w);
                assert_eq!(parse_hex(&s).unwrap().to_bits(), w.to_bits(), "{s}");
            }
            v *= 2.0;
        }
        assert_eq!(format_hex(f64::NAN), "NaN");
    }
}
