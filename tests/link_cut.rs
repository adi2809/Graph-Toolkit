use graph_library::structures::LinkCutForest;

#[test]
fn interleaved_operations_match_explicit_parent_paths() {
    let n = 35;
    let mut forest = LinkCutForest::new(n);
    let mut parent = vec![None::<(usize, u64)>; n];
    let mut state = 91u64;
    let mut rng = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        state
    };
    for i in 0..60_000 {
        let u = rng() as usize % n;
        let v = rng() as usize % n;
        let cap = if i % 29 == 0 { u64::MAX } else { rng() % 19 };
        let mut root = u;
        let mut path = Vec::new();
        while let Some((p, c)) = parent[root] {
            path.push((c, root));
            root = p;
        }
        match rng() % 5 {
            0 => {
                let mut top = v;
                while let Some((p, _)) = parent[top] {
                    top = p;
                }
                let can = parent[u].is_none() && top != u;
                assert_eq!(forest.link(u, v, cap), can);
                if can {
                    parent[u] = Some((v, cap));
                }
            }
            1 => {
                assert_eq!(forest.cut(u), parent[u].take().map(|x| x.1));
            }
            2 => {
                let ok = path.iter().all(|&(c, _)| c >= cap);
                assert_eq!(forest.subtract_path(u, cap), ok);
                if ok {
                    for (_, v) in path {
                        parent[v].as_mut().unwrap().1 -= cap;
                    }
                }
            }
            3 => assert_eq!(forest.root(u), root),
            _ => assert_eq!(forest.path_minimum(u), path.iter().min().copied()),
        }
    }
}
