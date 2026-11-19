//! Implicit rectangular grids with bit-packed obstacles and no allocated edges.

use std::error::Error;
use std::fmt;

use crate::shortest_path::astar_core;
use crate::traversal::{BreadthFirst, bfs_core};

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Cell {
    pub row: usize,
    pub col: usize,
}

impl Cell {
    pub const fn new(row: usize, col: usize) -> Self {
        Self { row, col }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Movement {
    Four,
    /// Unit-cost diagonals are enabled when both adjacent cardinal cells are open.
    Eight,
}

#[derive(Clone, Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum GridError {
    InvalidDimensions,
    OutOfBounds(Cell),
    Blocked(Cell),
}

impl fmt::Display for GridError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            Self::InvalidDimensions => f.write_str("grid dimensions must be nonzero and fit usize"),
            Self::OutOfBounds(cell) => {
                write!(f, "cell ({}, {}) is outside the grid", cell.row, cell.col)
            }
            Self::Blocked(cell) => write!(f, "cell ({}, {}) is blocked", cell.row, cell.col),
        }
    }
}

impl Error for GridError {}

#[derive(Clone, Debug)]
pub struct Grid {
    rows: usize,
    cols: usize,
    movement: Movement,
    blocked: Vec<u64>,
}

impl Grid {
    pub fn new(rows: usize, cols: usize, movement: Movement) -> Result<Self, GridError> {
        let size = rows
            .checked_mul(cols)
            .filter(|&size| size != 0)
            .ok_or(GridError::InvalidDimensions)?;
        Ok(Self {
            rows,
            cols,
            movement,
            blocked: vec![0; size.div_ceil(64)],
        })
    }

    pub fn dimensions(&self) -> (usize, usize) {
        (self.rows, self.cols)
    }

    pub fn set_blocked(&mut self, cell: Cell, blocked: bool) -> Result<(), GridError> {
        let index = self.index(cell).ok_or(GridError::OutOfBounds(cell))?;
        let mask = 1 << (index % 64);
        if blocked {
            self.blocked[index / 64] |= mask;
        } else {
            self.blocked[index / 64] &= !mask;
        }
        Ok(())
    }

    pub fn is_blocked(&self, cell: Cell) -> Result<bool, GridError> {
        let index = self.index(cell).ok_or(GridError::OutOfBounds(cell))?;
        Ok(self.blocked_at(index))
    }

    /// Returns distances, nearest sources, and paths in O(rows × cols) time.
    pub fn multi_source_bfs(&self, sources: &[Cell]) -> Result<GridDistances, GridError> {
        let sources = sources
            .iter()
            .map(|&cell| self.open_index(cell))
            .collect::<Result<Vec<_>, _>>()?;
        let forest = bfs_core(self.rows * self.cols, &sources, |index| {
            self.neighbors(index)
        });
        Ok(GridDistances {
            rows: self.rows,
            cols: self.cols,
            forest,
        })
    }

    /// Returns a shortest walk through open cells, including both endpoints.
    /// Uses Manhattan distance for cardinal movement and Chebyshev distance for diagonals.
    pub fn shortest_path(&self, from: Cell, to: Cell) -> Result<Option<Vec<Cell>>, GridError> {
        let source = self.open_index(from)?;
        let target = self.open_index(to)?;
        let route = astar_core(
            self.rows * self.cols,
            source,
            target,
            |index| self.neighbors(index).map(|next| (next, 0, 1)),
            |index| {
                let cell = cell_at(index, self.cols);
                let row = cell.row.abs_diff(to.row);
                let col = cell.col.abs_diff(to.col);
                match self.movement {
                    Movement::Four => (row + col) as u128,
                    Movement::Eight => row.max(col) as u128,
                }
            },
        );
        Ok(route.map(|path| {
            path.nodes
                .into_iter()
                .map(|id| cell_at(id, self.cols))
                .collect()
        }))
    }

    fn neighbors(&self, index: usize) -> impl Iterator<Item = usize> + '_ {
        const OFFSETS: [(isize, isize); 8] = [
            (-1, 0),
            (0, -1),
            (0, 1),
            (1, 0),
            (-1, -1),
            (-1, 1),
            (1, -1),
            (1, 1),
        ];
        let from = cell_at(index, self.cols);
        let count = if self.movement == Movement::Four {
            4
        } else {
            8
        };
        OFFSETS[..count].iter().filter_map(move |&(dr, dc)| {
            let row = from.row.checked_add_signed(dr)?;
            let col = from.col.checked_add_signed(dc)?;
            let next = self.index(Cell::new(row, col))?;
            if self.blocked_at(next)
                || (dr != 0
                    && dc != 0
                    && (self.blocked_at(row * self.cols + from.col)
                        || self.blocked_at(from.row * self.cols + col)))
            {
                return None;
            }
            Some(next)
        })
    }

    fn index(&self, cell: Cell) -> Option<usize> {
        index_at(cell, self.rows, self.cols)
    }

    fn open_index(&self, cell: Cell) -> Result<usize, GridError> {
        let index = self.index(cell).ok_or(GridError::OutOfBounds(cell))?;
        if self.blocked_at(index) {
            return Err(GridError::Blocked(cell));
        }
        Ok(index)
    }

    fn blocked_at(&self, index: usize) -> bool {
        self.blocked[index / 64] & (1 << (index % 64)) != 0
    }
}

#[derive(Clone, Debug)]
pub struct GridDistances {
    rows: usize,
    cols: usize,
    forest: BreadthFirst,
}

impl GridDistances {
    pub fn order(&self) -> impl ExactSizeIterator<Item = Cell> + '_ {
        self.forest
            .order()
            .iter()
            .map(|&index| cell_at(index, self.cols))
    }

    pub fn distance(&self, cell: Cell) -> Option<usize> {
        self.forest.distance(index_at(cell, self.rows, self.cols)?)
    }

    pub fn nearest_source(&self, cell: Cell) -> Option<Cell> {
        let source = self.forest.source(index_at(cell, self.rows, self.cols)?)?;
        Some(cell_at(source, self.cols))
    }

    pub fn path_to(&self, cell: Cell) -> Option<Vec<Cell>> {
        let path = self.forest.path_to(index_at(cell, self.rows, self.cols)?)?;
        Some(
            path.into_iter()
                .map(|index| cell_at(index, self.cols))
                .collect(),
        )
    }
}

fn cell_at(index: usize, cols: usize) -> Cell {
    Cell::new(index / cols, index % cols)
}

fn index_at(cell: Cell, rows: usize, cols: usize) -> Option<usize> {
    (cell.row < rows && cell.col < cols).then(|| cell.row * cols + cell.col)
}
