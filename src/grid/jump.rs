//! Horizontal-first jump point search for four-connected grids.
use super::{Cell, Grid, GridError, Movement, cell_at};
use crate::shortest_path::astar_core;

// Horizontal arrivals permit either vertical turn. A vertical arrival
// permits continuing in its direction and taking a forced horizontal turn.
const HORIZONTAL: usize = 0;
const NORTH: usize = 1;
const SOUTH: usize = 2;
const STATES: usize = 3;

impl Grid {
    fn open_step(&self, cell: Cell, dr: isize, dc: isize) -> Option<Cell> {
        let next = Cell::new(
            cell.row.checked_add_signed(dr)?,
            cell.col.checked_add_signed(dc)?,
        );
        (!self.blocked_at(self.index(next)?)).then_some(next)
    }

    fn jump_successor(
        &self,
        cell: Cell,
        arrival: usize,
        dr: isize,
        dc: isize,
        goal: Cell,
    ) -> Option<(usize, usize, u64)> {
        if arrival != HORIZONTAL {
            let incoming = if arrival == NORTH { -1 } else { 1 };
            if dr == -incoming {
                return None;
            }
            let previous = self.open_step(cell, -incoming, 0)?;
            if dc != 0 && self.open_step(previous, 0, dc).is_some() {
                return None;
            }
        }
        let mut current = cell;
        let mut length = 0;
        loop {
            let next = self.open_step(current, dr, dc)?;
            length += 1;
            let forced = dr != 0
                && [-1, 1].into_iter().any(|side| {
                    self.open_step(current, 0, side).is_none()
                        && self.open_step(next, 0, side).is_some()
                });
            if dc != 0 || next == goal || forced {
                let direction = if dc != 0 || next == goal {
                    HORIZONTAL
                } else if dr < 0 {
                    NORTH
                } else {
                    SOUTH
                };
                return Some((self.index(next)? * STATES + direction, 0, length));
            }
            current = next;
        }
    }

    /// Exact cardinal-grid path using JPS4. Vertical rays skip cells until
    /// the goal or an obstacle forces a turn. Returns each cell of the path.
    /// A grid configured for diagonals returns RequiresFourWayMovement.
    pub fn jump_point_path(&self, from: Cell, to: Cell) -> Result<Option<Vec<Cell>>, GridError> {
        if self.movement != Movement::Four {
            return Err(GridError::RequiresFourWayMovement);
        }
        let start = self.open_index(from)?;
        let goal = self.open_index(to)?;
        let states = (self.rows * self.cols)
            .checked_mul(STATES)
            .ok_or(GridError::InvalidDimensions)?;
        let route = astar_core(
            states,
            start * STATES,
            goal * STATES,
            |state| {
                let cell = cell_at(state / STATES, self.cols);
                [(0, -1), (0, 1), (-1, 0), (1, 0)]
                    .into_iter()
                    .filter_map(move |(dr, dc)| {
                        self.jump_successor(cell, state % STATES, dr, dc, to)
                    })
            },
            |state| {
                let cell = cell_at(state / STATES, self.cols);
                (cell.row.abs_diff(to.row) + cell.col.abs_diff(to.col)) as u128
            },
        );
        Ok(route.map(|path| {
            let mut cells = vec![from];
            let mut current = from;
            for state in path.nodes.into_iter().skip(1) {
                let end = cell_at(state / STATES, self.cols);
                while current != end {
                    if end.row > current.row {
                        current.row += 1;
                    } else if end.row < current.row {
                        current.row -= 1;
                    } else {
                        current.col = end.col;
                    }
                    cells.push(current);
                }
            }
            cells
        }))
    }
}
