use std::path::PathBuf;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use super::Kmer;
use std::collections::HashMap;
use petgraph::graph::NodeIndex;

// this module takes care of outputting the necessary information
// for reconstructing sequences in assemblies. otherwise having just 
// the graph would be a bit useless

pub fn write_minimizers_and_seq_of_kmers(output_prefix :&PathBuf, node_indices :&HashMap<Kmer,NodeIndex>, kmer_seqs :&HashMap<Kmer,String>, kmer_origin: &HashMap<Kmer,String>, dbg_nodes :&HashMap<Kmer,u32> , k:usize, l:usize)
{
    let output_filename = format!("{}{}",output_prefix.to_string_lossy(),".sequences");
    let path = PathBuf::from(&output_filename);
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };

    write!(file, "# k = {}\n",k).unwrap();
    write!(file, "# l = {}\n",l).unwrap();
    write!(file, "# structure of remaining of the file:\n").unwrap();
    write!(file, "# [node name]\t[list of minimizers]\t[sequence of node]\t[abundance]\t[origin]\n").expect("error writing sequences file");

    for (node, index) in node_indices {
        let seq       = &kmer_seqs[node];
        let abundance = &dbg_nodes[node];
        let origin  = &kmer_origin[node];
        let s_line = format!("{}\t{}\t{}\t{}\t{}\n",index.index(),node.print_as_string(), seq, abundance, origin);
        write!(file, "{}", s_line).expect("error writing minimizers/sequences");
    }


}

