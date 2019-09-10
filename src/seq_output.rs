use std::path::Path;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use super::Kmer;
use std::collections::HashMap;
use petgraph::graph::NodeIndex;

// this module takes care of outputting the necessary information
// for reconstructing sequences in assemblies. otherwise having just 
// the graph would be a bit useless

pub fn write_minimizers_and_seq_of_kmers(graph_filename :&str, node_indices :&HashMap<Kmer,NodeIndex>, kmer_seqs :&HashMap<Kmer,String>, k:usize, l:usize)
{
    let output_filename = format!("{}{}",Path::new(graph_filename).file_stem().unwrap().to_string_lossy(),".sequences");
    let path = Path::new(&output_filename);
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };

    write!(file, "# k = {}\n",k).unwrap();
    write!(file, "# l = {}\n",l).unwrap();
    write!(file, "# structure of remaining of the file:\n").unwrap();
    write!(file, "# [node name]\t[list of minimizers]\t[sequence of node]\n").expect("error writing sequences file");

    for (node, index) in node_indices {
        let seq = &kmer_seqs[node];
        let s_line = format!("{}\t{}\t{}\n",index.index(),node.print_as_string(),seq);
        write!(file, "{}", s_line).expect("error writing minimizers/sequences");
    }


}

