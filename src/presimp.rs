use petgraph::graph::DiGraph;
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use super::Kmer;
use crate::kmer_vec::get;
use std::collections::HashMap;

/// for each node v, examines its out neighbors
/// find the neighbor with maximum abundance M
/// consider also the abundance m of v
/// delete all edges to nodes having abundance less than
/// min(M,m)*some_factor where some_factor can be 0.1
pub fn presimp(gr: &mut DiGraph::<Kmer,Kmer>, dbg_nodes: &HashMap<Kmer,u32>, factor: f32)  {
    let nodes_vect : Vec<&Kmer> = dbg_nodes.keys().collect();

    for node in gr.node_indices() {
        let idx = node.index();
        let abundance = *dbg_nodes.get(nodes_vect[idx]).unwrap();

        //for dir in [Direction::Outgoing, Direction::Incoming].iter()
        //{
        let dir = &Direction::Outgoing;

            let mut abundances = Vec::new();
            let mut neighbors = gr.neighbors_directed(node, *dir).detach();
            while let Some(neigh) = neighbors.next_node(&gr) {
                let kmer2 = &nodes_vect[neigh.index()];
                let abundance2 = dbg_nodes[kmer2];
                abundances.push(abundance2);
            }
            if abundances.len() < 2 { continue; } // preserves connectivity
            let abundance_ref = std::cmp::min(*abundances.iter().max().unwrap(),abundance);

            let mut neighbors = gr.neighbors_directed(node, *dir).detach();
            while let Some(neigh) = neighbors.next_node(&gr) {
                let kmer2 = &nodes_vect[neigh.index()];
                let abundance2 = dbg_nodes[kmer2];
                if (abundance2 as f32) < factor*(abundance_ref as f32)
                {
                    let edge = if *dir == Direction::Incoming { gr.find_edge(neigh,node) } else { gr.find_edge(node, neigh) }; 
                    gr.remove_edge(edge.unwrap());
                }
            }
        //}

        // normalize so that if v->w was delete, delete w->v too
        let mut to_del = Vec::new();
        for e in gr.edge_references() {
            let edge = gr.find_edge(e.target(),e.source());
            if !edge.is_some()
            {
                to_del.push(e.id());
            }
        }
        for edge in to_del
        {
            gr.remove_edge(edge);
        }

    }
}
