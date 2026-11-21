//! Standalone comparison using identical seeded sparse and dense inputs.
#[path = "../tests/support/flow_certificate.rs"]
mod certificate;
#[path = "support/flow_instances.rs"]
mod inputs;
use graph_library::flow::{
    TreeCycleBasis, max_flow, max_flow_core_cycles, max_flow_dinic, max_flow_tree_cycle,
};
use std::time::Instant;

fn main() {
    const SOLVERS: [&str; 4] = [
        "dinic",
        "push_relabel",
        "core_cycles",
        "signed_barrier_core",
    ];
    println!(
        "family,n,arcs,solver,median_ms,min_ms,max_ms,value,epochs,augmentations,fractional_edges,barrier_steps,warm_flow"
    );
    for (family, n) in [
        ("sparse", 16384),
        ("sparse", 65536),
        ("dense", 512),
        ("dense", 1024),
    ] {
        let g = inputs::instance(family, n);
        let mut times = vec![Vec::new(); SOLVERS.len()];
        let mut details = vec![None; SOLVERS.len()];
        let mut optimum = None;
        for round in 0..6 {
            for offset in 0..SOLVERS.len() {
                let i = (round + offset) % SOLVERS.len();
                let start = Instant::now();
                let mut stats = [None; 5];
                let flow = match i {
                    0 => max_flow_dinic(&g, 0, n - 1).unwrap(),
                    1 => max_flow(&g, 0, n - 1).unwrap(),
                    2 => {
                        let r = max_flow_core_cycles(&g, 0, n - 1).unwrap();
                        stats[..3].copy_from_slice(
                            &[r.epochs, r.augmentations, r.fractional_edges]
                                .map(|v| Some(v as u128)),
                        );
                        r.flow
                    }
                    _ => {
                        let r =
                            max_flow_tree_cycle(&g, 0, n - 1, 8, TreeCycleBasis::DynamicHierarchy)
                                .unwrap();
                        stats[3] = Some(r.descent_steps as u128);
                        stats[4] = Some(r.value_before_completion);
                        r.flow
                    }
                };
                let ms = start.elapsed().as_secs_f64() * 1000.;
                certificate::verify(&g, 0, n - 1, &flow);
                assert_eq!(*optimum.get_or_insert(flow.value), flow.value);
                assert_eq!(*details[i].get_or_insert(stats), stats);
                if round > 0 {
                    times[i].push(ms);
                }
            }
        }
        for (i, mut sample) in times.into_iter().enumerate() {
            sample.sort_by(f64::total_cmp);
            let stats = details[i]
                .unwrap()
                .map(|v| v.map_or_else(|| "-".into(), |x| x.to_string()))
                .join(",");
            println!(
                "{family},{n},{},{},{:.3},{:.3},{:.3},{},{stats}",
                g.edge_count(),
                SOLVERS[i],
                sample[2],
                sample[0],
                sample[4],
                optimum.unwrap()
            );
        }
    }
}
