use graph_library::grid::{Cell, Grid, GridError, Movement};
#[test]
fn all_three_by_three_maps_and_open_terminal_pairs() {
    for mask in 0..512 {
        let mut g = Grid::new(3, 3, Movement::Four).unwrap();
        let cell = |i: usize| Cell::new(i / 3, i % 3);
        for i in 0..9 {
            g.set_blocked(cell(i), mask >> i & 1 != 0).unwrap();
        }
        for i in 0..9 {
            if mask >> i & 1 != 0 {
                continue;
            }
            let bfs = g.multi_source_bfs(&[cell(i)]).unwrap();
            for j in 0..9 {
                if mask >> j & 1 != 0 {
                    continue;
                }
                let path = g.jump_point_path(cell(i), cell(j)).unwrap();
                assert_eq!(path.as_ref().map(|p| p.len() - 1), bfs.distance(cell(j)));
                if let Some(p) = path {
                    assert_eq!((p.first(), p.last()), (Some(&cell(i)), Some(&cell(j))));
                    assert!(p.iter().all(|&c| !g.is_blocked(c).unwrap()));
                    for a in p.windows(2) {
                        assert_eq!(a[0].row.abs_diff(a[1].row) + a[0].col.abs_diff(a[1].col), 1);
                    }
                }
            }
        }
    }
    let g = Grid::new(1, 1, Movement::Eight).unwrap();
    assert_eq!(
        g.jump_point_path(Cell::new(0, 0), Cell::new(0, 0)),
        Err(GridError::RequiresFourWayMovement)
    );
}
