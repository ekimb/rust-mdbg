use petgraph::graph::DiGraph;
use petgraph::visit::EdgeRef;
use std::path::Path;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use super::Kmer;
use std::collections::HashMap;

fn determine_orientation(id1: usize, id2: usize, nodes_vect: &Vec<&Kmer> ) -> (&'static str,&'static str)
{
    let n1 = nodes_vect[id1];
    let n2 = nodes_vect[id2];
    let rev_n1 = n1.reverse();
    let rev_n2 = n2.reverse();
    if n1.prefix() == n2.suffix() { return ("+","+"); }
    if n1.prefix() == rev_n2.suffix() { return ("+","-"); }
    if rev_n1.prefix() == n2.suffix() { return ("-","+"); }
    if rev_n1.prefix() == rev_n2.suffix() { return ("-","-"); }
    panic!("unknown orientation");
}   

pub fn output_gfa(gr: &DiGraph::<Kmer,Kmer>, dbg_nodes: &HashMap<Kmer,u32>) {
    let nodes_vect : Vec<&Kmer> = dbg_nodes.keys().collect();
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
