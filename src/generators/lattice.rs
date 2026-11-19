use super::{add, null_graph};
use crate::{Graph, GraphBuilder, GraphError, GraphKind, NONE};

/// A triangular lattice of `rows` by `columns` triangles. Vertex IDs are
/// row-major. Periodic dimensions require at least three rows and five columns.
pub fn triangular_lattice_graph(
    rows: usize,
    columns: usize,
    periodic: bool,
) -> Result<Graph<()>, GraphError> {
    if rows == 0 || columns == 0 {
        return Ok(null_graph());
    }
    if periodic && (rows < 3 || columns < 5) {
        return Err(GraphError::InvalidModelParameters);
    }
    let width = columns.div_ceil(2) + usize::from(!periodic);
    let height = rows
        .checked_add(usize::from(!periodic))
        .ok_or(GraphError::GraphSizeOverflow)?;
    let omit = !periodic && columns % 2 != 0;
    let n = width
        .checked_mul(height)
        .ok_or(GraphError::GraphSizeOverflow)?
        - if omit { height / 2 } else { 0 };
    let row_width = |y| width - usize::from(omit && y % 2 != 0);
    let node = |x, y| y * width + x - if omit { y / 2 } else { 0 };
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    for y in 0..height {
        let upper = (y + 1) % height;
        for x in 0..row_width(y) {
            let from = node(x, y);
            if x + 1 < row_width(y) || periodic {
                add(&mut builder, from, node((x + 1) % width, y));
            }
            if upper == 0 && !periodic {
                continue;
            }
            if x < row_width(upper) {
                add(&mut builder, from, node(x, upper));
            }
            let diagonal = if y % 2 != 0 {
                (x + 1 < row_width(upper) || periodic).then_some((x + 1) % width)
            } else if x == 0 {
                periodic.then_some(width - 1)
            } else {
                Some(x - 1)
            };
            if let Some(x) = diagonal {
                add(&mut builder, from, node(x, upper));
            }
        }
    }
    Ok(builder.build())
}

/// A honeycomb lattice of `rows` by `columns` hexagons. Vertex IDs are
/// row-major. Periodic lattices require at least two rows and an even number
/// of columns greater than one.
pub fn hexagonal_lattice_graph(
    rows: usize,
    columns: usize,
    periodic: bool,
) -> Result<Graph<()>, GraphError> {
    if rows == 0 || columns == 0 {
        return Ok(null_graph());
    }
    if periodic && (rows < 2 || columns < 2 || columns % 2 != 0) {
        return Err(GraphError::InvalidModelParameters);
    }
    let width = columns
        .checked_add(usize::from(!periodic))
        .ok_or(GraphError::GraphSizeOverflow)?;
    let height = rows
        .checked_mul(2)
        .and_then(|n| n.checked_add(if periodic { 0 } else { 2 }))
        .ok_or(GraphError::GraphSizeOverflow)?;
    let count = width
        .checked_mul(height)
        .ok_or(GraphError::GraphSizeOverflow)?;
    let mut ids = Vec::with_capacity(count);
    let mut n = 0;
    for y in 0..height {
        for x in 0..width {
            if !periodic
                && ((x == 0 && y == height - 1)
                    || (x == columns && y == if columns % 2 == 0 { 0 } else { height - 1 }))
            {
                ids.push(NONE);
            } else {
                ids.push(n);
                n += 1;
            }
        }
    }
    let mut builder = GraphBuilder::new(n, GraphKind::Undirected);
    for y in 0..height {
        for x in 0..width {
            let from = ids[y * width + x];
            if from == NONE {
                continue;
            }
            if y + 1 < height || periodic {
                let to = ids[((y + 1) % height) * width + x];
                if to != NONE {
                    add(&mut builder, from, to);
                }
            }
            if x % 2 == y % 2 && (x + 1 < width || periodic) {
                let to = ids[y * width + (x + 1) % width];
                if to != NONE {
                    add(&mut builder, from, to);
                }
            }
        }
    }
    Ok(builder.build())
}
