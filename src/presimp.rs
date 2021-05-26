use petgraph::graph::DiGraph;
use petgraph::graph::{EdgeIndex};
use petgraph::Direction;
use super::Kmer;
use crate::kmer_vec::get;
use std::collections::HashMap;
use super::thread_update_vec;
use std::sync::{Arc, Mutex, MutexGuard};
use crossbeam_utils::{thread};
use closure::closure;

/// For each node v, examine outgoing neighbors;
/// Find the neighbor with maximum abundance M;
/// Consider also the abundance m of v;
/// Delete all edges to nodes having abundance less than min(M,m)*f, e.g., f = 0.1.

pub fn find_removed_edges(gr: &DiGraph::<Kmer, Kmer>, dbg_nodes: &HashMap<Kmer, u32>, factor: f32, threads: usize) -> Vec<EdgeIndex> {
    let mut chunks = gr.node_indices().collect::<Vec<_>>();
    let mut chunk_length = 1;
    let mut removed_edges_all = Arc::new(Mutex::new(HashMap::<usize, Vec<EdgeIndex>>::new()));
    let nodes_vect : Vec<&Kmer> = dbg_nodes.keys().collect();
    if chunks.len() > threads {chunk_length = chunks.len()/threads+1;}
    thread::scope(|s| {
        let mut guards = Vec::with_capacity(threads);
        for (thread_num, chunk) in chunks.chunks(chunk_length).enumerate() {
            let removed_edges_all = removed_edges_all.clone();
            let guard = s.spawn(closure!(move chunk, ref gr, ref dbg_nodes, ref factor, ref nodes_vect, |_| {
                let mut removed_edges = Vec::<EdgeIndex>::new();
                for node in chunk {
                    let idx = node.index();
                    let abundance = *dbg_nodes.get(nodes_vect[idx]).unwrap();
                    let dir = &Direction::Outgoing;
                    let mut abundances = Vec::new();
                    let mut neighbors = gr.neighbors_directed(*node, *dir).detach();
                    while let Some(neigh) = neighbors.next_node(&gr) {
                        let kmer2 = &nodes_vect[neigh.index()];
                        let abundance2 = dbg_nodes[kmer2];
                        abundances.push(abundance2);
                    }
                    if abundances.len() < 2 {continue;} // preserves connectivity
                    let abundance_ref = std::cmp::min(*abundances.iter().max().unwrap(),abundance);
                    let mut neighbors = gr.neighbors_directed(*node, *dir).detach();
                    while let Some(neigh) = neighbors.next_node(&gr) {
                        let kmer2 = &nodes_vect[neigh.index()];
                        let abundance2 = dbg_nodes[kmer2];
                        if (abundance2 as f32) < factor * (abundance_ref as f32) {
                            let edge = if *dir == Direction::Incoming {gr.find_edge(neigh, *node).unwrap()} else {gr.find_edge(*node, neigh).unwrap()}; 
                            removed_edges.push(edge);
                            let other_edge = gr.find_edge(gr.raw_edges()[edge.index()].target(), gr.raw_edges()[edge.index()].source());
                            if other_edge.is_some() {removed_edges.push(other_edge.unwrap())}
                        }
                    }            
                
                }
                thread_update_vec(&removed_edges_all, removed_edges, thread_num);
            }));
            guards.push(guard);
        }
    }).unwrap();
    let mut edges = Vec::<EdgeIndex>::new();
    let mut removed_edges_all = removed_edges_all.lock().unwrap();
    for thread_num in 0..threads {
        let mut removed_edges = removed_edges_all.entry(thread_num).or_insert(Vec::new());
        for edge in removed_edges.iter() {
            edges.push(edge.clone())
        }
    }
    edges
}

pub fn presimp(gr: &mut DiGraph::<Kmer, Kmer>, edges: &Vec<EdgeIndex>) {
    for edge in edges.iter() {
        gr.remove_edge(*edge);
    }
}