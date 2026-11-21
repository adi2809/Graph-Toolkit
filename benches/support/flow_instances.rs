use graph_library::{Graph, GraphKind};
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn capacity(&mut self) -> u64 {
        1 + self.next() % 1000
    }
}
pub fn instance(family: &str, n: usize) -> Graph<u64> {
    let mut r = Rng(0x73d4bbe9);
    let mut g = Graph::builder(n, GraphKind::Directed);
    match family {
        "sparse" => {
            for u in 0..n - 1 {
                g.add_edge(u, u + 1, r.capacity()).unwrap();
            }
            for _ in 0..5 * n {
                let u = r.next() as usize % n;
                let v = r.next() as usize % n;
                g.add_edge(u, v, r.capacity()).unwrap();
            }
            for v in (1..n - 1).step_by(256) {
                g.add_edge(0, v, r.capacity()).unwrap();
                g.add_edge(v, n - 1, r.capacity()).unwrap();
            }
        }
        "dense" => {
            for u in 0..n {
                for v in 0..n {
                    if u != v && r.next() % 2 == 0 {
                        g.add_edge(u, v, r.capacity()).unwrap();
                    }
                }
            }
        }
        "grid" => {
            let side = n.isqrt();
            assert_eq!(side * side, n);
            for u in 0..n {
                for v in [
                    (u % side + 1 < side).then_some(u + 1),
                    (u + side < n).then_some(u + side),
                ]
                .into_iter()
                .flatten()
                {
                    g.add_edge(u, v, r.capacity()).unwrap();
                    g.add_edge(v, u, r.capacity()).unwrap();
                }
            }
        }
        _ => unreachable!(),
    }
    g.build()
}
