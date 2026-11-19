#![forbid(unsafe_code)]
#![doc = include_str!("../README.md")]

mod graph;
mod queue;
mod random;

pub mod centrality;
pub mod clustering;
pub mod community;
pub mod connectivity;
pub mod flow;
pub mod functions;
pub mod generators;
pub mod grid;
pub mod isomorphism;
pub mod matching;
#[cfg(feature = "multicommodity")]
pub mod multicommodity;
pub mod planarity;
pub mod shortest_path;
pub mod traversal;

pub use graph::{Arc, Edge, EdgeId, Graph, GraphBuilder, GraphError, GraphKind, NodeId};

pub(crate) const NONE: usize = usize::MAX;
