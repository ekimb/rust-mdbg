use std::path::PathBuf;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use itertools::Itertools;
use std::io::{BufWriter, BufRead, BufReader};
use super::Kmer;
use std::collections::HashMap;
use petgraph::graph::NodeIndex;
use std::io::{self};
use std::path::Path;


// this module takes care of outputting the necessary information
// for reconstructing sequences in assemblies. otherwise having just 
// the graph would be a bit useless

pub fn write_minimizers_and_seq_of_kmers(output_prefix :&PathBuf, node_indices :&mut HashMap<Kmer,usize>, kmer_origin: &HashMap<Kmer,String>, dbg_nodes :&HashMap<Kmer,u32> , k:usize, l:usize)
{
    let output_filename = format!("{}{}",output_prefix.to_string_lossy(),".sequences");
    let path = PathBuf::from(&output_filename);
        let mut out_file = match File::create(path) {
            Err(why) => panic!("couldn't load sequences file: {}", why.description()),
            Ok(out_file) => out_file,
    }; 
    write!(out_file, "# k = {}\n",k).unwrap();
    write!(out_file, "# l = {}\n",l).unwrap();
    write!(out_file, "# structure of remaining of the file:\n").unwrap();
    write!(out_file, "# [node name]\t[list of minimizers]\t[sequence of node]\t[abundance]\t[origin]\t[shift]\n").expect("error writing sequences file");
    let mut node_count = 0;
    let input_filename = format!("{}.0.sequences",output_prefix.to_string_lossy());
    let input_path = PathBuf::from(&input_filename);
    if let Ok(lines) = read_lines(input_path) {
        for line in lines {
            if let Ok(ip) = line {
                //println!("{}", ip);
                let str_vec : Vec<&str> = ip.trim().split("\t").collect();
                let node_raw : Vec<u64> = str_vec[0].split(|c| c == '[' || c == ']' || c == ' ' || c == ',').filter(|&c| c != "").map(|s| s.parse().unwrap()).collect::<Vec<u64>>();
                let kmer_seq = str_vec[1].to_string();  
                let kmer_origin = str_vec[2].to_string(); 
                let node_shift : Vec<usize> = str_vec[3].trim().split(|c| c == '(' || c == ')' || c == ' ' || c == ',').filter(|&c| c != "").map(|s| s.parse().unwrap()).collect::<Vec<usize>>();
                let mut node : Kmer = Kmer::make_from(&node_raw);
                let max_node = usize::max_value();
                let index = node_indices.get(&node).unwrap_or(&max_node);
                if index != &max_node {
                    let abundance = dbg_nodes[&node];
                    node_count += 1;
                    let s_line = format!("{}\t{}\t{}\t{}\t{}\t{:?}",index,node.print_as_string(), kmer_seq, abundance, kmer_origin, node_shift);
                    write!(out_file, "{}\n", s_line).expect("error writing minimizers/sequences");
                    node_indices.remove(&node);
                }

            }
        }
    }
   /* for (node, index) in node_indices {
        println!("{} {}", node.print_as_string(), dbg_nodes[&node])

    }*/

}
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

