use petgraph::graph::DiGraph;
use petgraph::visit::EdgeRef;
use std::path::Path;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use super::k;
use super::reverse_node;
use std::collections::HashMap;

fn determine_orientation(id1: usize, id2: usize, nodes_vect: &Vec<&[u32;k]> ) -> (&'static str,&'static str)
{
    let n1 = nodes_vect[id1];
    let n2 = nodes_vect[id2];
    let rev_n1 = reverse_node(*n1);
    let rev_n2 = reverse_node(*n2);
    if n1[..k-1] == n2[1..k] { return ("+","+"); }
    if n1[..k-1] == rev_n2[1..k] { return ("+","-"); }
    if rev_n1[..k-1] == n2[1..k] { return ("-","+"); }
    if rev_n1[..k-1] == rev_n2[1..k] { return ("-","-"); }
    panic!("unknown orientation");
}   

pub fn output_gfa(gr: &DiGraph::<[u32;k],[u32;k]>, dbg_nodes: &HashMap<[u32;k],u32>) {
    let nodes_vect : Vec<&[u32;k]> = dbg_nodes.keys().collect();
    let path = Path::new("graph.gfa");
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };
    write!(file, "H\tVZ:Z:1\n").expect("error writing GFA header");

    for node in gr.node_indices() {
        let s_line = format!("S\t{}\t*\tLN:i:100\n",node.index());
        write!(file, "{}", s_line).expect("error writing s_line");
    }

    for e in gr.edge_references() {
        let id1 = e.source().index();
        let id2 = e.target().index();
        let (ori1, ori2) = determine_orientation(id1,id2,&nodes_vect);
        let l_line = format!("L\t{}\t{}\t{}\t{}\t50M\n", id1, ori1, id2, ori2);
        write!(file, "{}", l_line).expect("error writing l_line");
    }
}
