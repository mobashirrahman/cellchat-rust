//! CellChat's ComputeSNN: sparse shared-neighbour counts and Jaccard pruning.
use std::collections::BTreeMap;

#[derive(Debug)]
pub struct SparseSnn {
    pub i: Vec<i32>,
    pub p: Vec<i32>,
    pub x: Vec<f64>,
}

/// `nn_ranked` is a column-major matrix of one-based neighbour indices.
/// Duplicate neighbours contribute multiplicities, as Eigen's setFromTriplets does.
pub fn compute_snn(
    nn_ranked: &[i32],
    rows: usize,
    cols: usize,
    prune: f64,
) -> Result<SparseSnn, String> {
    if nn_ranked.len() != rows * cols {
        return Err("neighbour matrix dimensions do not match its contents".into());
    }
    let mut memberships = vec![BTreeMap::<usize, f64>::new(); rows];
    for col in 0..cols {
        for row in 0..rows {
            let neighbour = nn_ranked[col * rows + row];
            if neighbour < 1 || neighbour as usize > rows {
                return Err("neighbour index is out of bounds".into());
            }
            *memberships[neighbour as usize - 1].entry(row).or_default() += 1.0;
        }
    }
    let mut overlaps = vec![BTreeMap::<usize, f64>::new(); rows];
    for members in memberships {
        for (&col, &right) in &members {
            for (&row, &left) in &members {
                *overlaps[col].entry(row).or_default() += left * right;
            }
        }
    }
    let mut out = SparseSnn {
        i: Vec::new(),
        p: vec![0],
        x: Vec::new(),
    };
    let k = cols as f64;
    for column in overlaps {
        for (row, overlap) in column {
            let value = overlap / (k + (k - overlap));
            let value = if value < prune { 0.0 } else { value };
            // Eigen removes explicit zero entries, retaining non-finite values.
            if value != 0.0 {
                out.i.push(row as i32);
                out.x.push(value);
            }
        }
        out.p.push(out.i.len() as i32);
    }
    Ok(out)
}
