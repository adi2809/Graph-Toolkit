//! Fixed seeded instances, rotated run order, certificates outside the timer.
#[path = "../tests/support/flow_certificate.rs"]
mod certificate;
use graph_library::flow::{
    TreeCycleBasis, max_flow, max_flow_dinic, max_flow_dynamic_trees, max_flow_st_planar,
    max_flow_tree_cycle,
};
use std::time::Instant;

#[path = "support/flow_instances.rs"]
mod flow_instances;
use flow_instances::instance;
fn main() {
    println!("family,n,m,method,median_ms,min_ms,max_ms,value,warm_flow,steps");
    for (family, n) in [
        ("sparse", 16_384),
        ("sparse", 65_536),
        ("dense", 512),
        ("dense", 1024),
        ("grid", 128 * 128),
    ] {
        let g = instance(family, n);
        let count = if family == "grid" { 6 } else { 5 };
        let mut times = vec![Vec::new(); count];
        let mut expected = None;
        let mut warm = [(0, 0); 2];
        // One warm-up and five measured runs per solver.
        for round in 0..6 {
            for k in 0..count {
                let i = (k + round) % count;
                let start = Instant::now();
                let flow = match i {
                    0 => max_flow_dinic(&g, 0, n - 1).unwrap(),
                    1 => max_flow(&g, 0, n - 1).unwrap(),
                    2 | 3 => {
                        let basis = [
                            TreeCycleBasis::ResidualBfs,
                            TreeCycleBasis::ResidualAndLowStretch,
                        ][i - 2];
                        let f = max_flow_tree_cycle(&g, 0, n - 1, 8, basis).unwrap();
                        warm[i - 2] = (f.value_before_completion, f.descent_steps);
                        f.flow
                    }
                    4 => max_flow_dynamic_trees(&g, 0, n - 1).unwrap(),
                    _ => max_flow_st_planar(&g, 0, n - 1).unwrap().unwrap(),
                };
                let elapsed = start.elapsed().as_secs_f64() * 1000.;
                certificate::verify(&g, 0, n - 1, &flow);
                assert_eq!(*expected.get_or_insert(flow.value), flow.value);
                if round > 0 {
                    times[i].push(elapsed);
                }
            }
        }
        for (i, mut sample) in times.into_iter().enumerate() {
            sample.sort_by(f64::total_cmp);
            let method = [
                "dinic",
                "push_relabel",
                "bfs_cycle_core",
                "low_stretch_cycle_core",
                "dynamic_trees",
                "planar_dual",
            ][i];
            let (warm_value, steps) = if i == 2 || i == 3 {
                (warm[i - 2].0.to_string(), warm[i - 2].1.to_string())
            } else {
                ("-".into(), "-".into())
            };
            println!(
                "{family},{n},{},{method},{:.3},{:.3},{:.3},{},{warm_value},{steps}",
                g.edge_count(),
                sample[2],
                sample[0],
                sample[4],
                expected.unwrap()
            );
        }
    }
}
