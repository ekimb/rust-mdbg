#![allow(unused_variables)]
#![allow(non_upper_case_globals)]

use pbr::ProgressBar;
use bio::io::fasta;
use std::collections::HashSet;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use fasthash::city;
use petgraph_graphml::GraphMl;
use petgraph::graph::DiGraph;
use petgraph::graph::NodeIndex;
//use adler32::RollingAdler32;
use std::fs;
use structopt::StructOpt;
use std::path::PathBuf;

mod utils;
mod gfa_output;
mod seq_output;

mod kmer_vec;
// mod kmer_array; // not working yet

const revcomp_aware: bool = true; // shouldn't be set to false except for strand-directed data or for debugging

//use typenum::{U31,U32}; // for KmerArray
type Kmer = kmer_vec::KmerVec;
type Overlap= kmer_vec::KmerVec;

 
fn extract_minimizers(seq: &str, l :usize, percentage_retain_hashes :f64) -> Vec<String>
{
    minhash(seq, l, percentage_retain_hashes)
    //wk_minimizers(seq, percentage_retain_hashes) // unfinished
}

fn minhash(seq: &str, l:usize, percentage_retain_hashes :f64) -> Vec<String>
{
    //let hash_mode = 0; // 1: rolling (currently incompat with revcomp), 0: cityhash
    let mut res = Vec::new();
    if seq.len() < l { return res; }

    let mut nl =  4f32.powf(l as f32) as u32;
    if revcomp_aware { nl = nl / 2; }
    //let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
    let mut h0 :u32 = 0;
    for start in 0..(&seq.len()-l+1)
    {
        let mut lmer = String::from(&seq[start..start+l]);
        if revcomp_aware {
            let lmer_rev = utils::revcomp(&lmer);
            lmer = std::cmp::min(lmer, lmer_rev);
        }

        if true
        //if hash_mode == 0 
        {
            h0 = city::hash32(&lmer) % nl;
        }/*
        else if hash_mode == 1 {
            if start > 0
            {
                h1.remove(l, seq.as_bytes()[start-1]);
                h1.update(seq.as_bytes()[start + l-1]);
            }
        }*/
       
        //let considered_hash = if hash_mode == 0 { h0 }  else { h1.hash() % nl };
        let considered_hash = h0;
        if considered_hash < (nl as f64*percentage_retain_hashes) as u32 {
            res.push(lmer.to_string());
            //println!("selected lmer: {} (hash {})", lmer.to_string(), considered_hash);
        }
    }
    res
}

// TODO finish implementing it
#[allow(dead_code)]
fn wk_minimizers(seq: &str, l :usize) -> Vec<String>
{
    // actually it's (w,l)-minimizers, k is already used in the dbg
    let w = 15;
    let hash_mode = 0; // 0: cityhash
    let mut res = Vec::new();
    if seq.len() < l { return res; }
    let mut h0 :u32 = 0;
    let nl =  4f32.powf(l as f32) as u32;
    // TODO create a set of (minimizers,positions)
    
    for start in 0..(&seq.len()-w+1)
    {
        let mut min : u32 = std::u32::MAX;
        for pos in start..start+w
        {
            let mut lmer = String::from(&seq[start..start+l]);
            if revcomp_aware {
                let lmer_rev = utils::revcomp(&lmer);
                lmer = std::cmp::min(lmer, lmer_rev);
            }
            h0 = city::hash32(lmer) % nl;
            min = std::cmp::min(h0,min); // TODO record position as well
            // TODO add(minimizer,position) to set
        }

        //&res.push(lmer.to_string());
        res.push(h0.to_string()); // placeholder
        h0 = 1; // placeholder
    }

    // TODO convert the list of (minimizers,positions) to a vector of list minimizers
    res
}

#[derive(Debug, StructOpt)]
#[structopt(name = "rust-mhdbg", about = "Original implementation of Min-Hash de Bruijn graphs")]
struct Opt {
    /// Activate debug mode
    // short and long flags (-d, --debug) will be deduced from the field's name
    #[structopt(short, long)]
    debug: bool,

    /// Input file
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,

    #[structopt(short, long)]
    k: Option<usize>,
    #[structopt(short, long)]
    l: Option<usize>,
    #[structopt(short, long)]
    pch: Option<f64>,
    #[structopt(long)]
    test1: bool,
    #[structopt(long)]
    test2: bool,
}


fn main() {
    let opt = Opt::from_args();      

    let mut filename = PathBuf::new();
    let mut k: usize = 10;
    let mut l: usize = 12;
    let mut percentage_retain_hashes :f64 = 0.10;

    if !opt.reads.is_none() { filename = opt.reads.unwrap(); }
   
    if opt.test1 {
        filename = PathBuf::from("../read50x_ref10K_e001.fa"); 
        if opt.k.is_none() { k = 5; }
        if opt.l.is_none() { l = 8; }
        if opt.pch.is_none() { percentage_retain_hashes = 0.15 };
    }
    else if opt.test2 {
        filename = PathBuf::from("../SRR9969842_vs_chr4.fasta");
        if opt.k.is_none() { k = 50; }
        if opt.l.is_none() { l = 12; }
        if opt.pch.is_none() { percentage_retain_hashes = 0.01 };
    }
    
    if !opt.k.is_none() { k = opt.k.unwrap() } else { println!("Warning: using default k value ({})",k); } 
    if !opt.l.is_none() { l = opt.l.unwrap() } else { println!("Warning: using default l value ({})",l); }
    if !opt.pch.is_none() { percentage_retain_hashes = opt.pch.unwrap() } else { println!("Warning: using default hash rate ({}%)",percentage_retain_hashes*100.0); }

    let debug = opt.debug;

    // init some useful objects
    let mut somewhat_reads : Vec<Vec<String>> = vec![];
    let mut all_minimizers : HashSet<String> = HashSet::new();
    let mut nb_minimizers_per_read : f64 = 0.0;
    let mut nb_reads : u64 = 0;


    // get file size for progress bar
    let metadata = fs::metadata(&filename).expect("error opening input file");
    let file_size = metadata.len();
    let mut pb = ProgressBar::new(file_size);

    // fasta parsing
    // possibly swap it later for https://github.com/aseyboldt/fastq-rs
    let reader = fasta::Reader::from_file(&filename).unwrap();
    for result in reader.records() {
	    let record = result.unwrap();
	    let seq = record.seq();
	
        let seq_str = String::from_utf8_lossy(seq);

        //println!("seq: {}", seq_str);

        let minimizers = extract_minimizers(&seq_str, l, percentage_retain_hashes);

        somewhat_reads.push(minimizers.clone());
        for minim in minimizers.iter() {
            all_minimizers.insert(minim.to_string());
        }

        // stats
        nb_minimizers_per_read += minimizers.len() as f64;
        nb_reads += 1;
        pb.add(record.seq().len() as u64 + record.id().len() as u64); // get approx size of entry
        
    }
    pb.finish_print("done converting reads to minimizers");
    nb_minimizers_per_read /= nb_reads as f64;

    println!("avg number of minimizers/read: {}",nb_minimizers_per_read);

    // assign numbers to minimizers
    let mut minimizer_to_int : HashMap<String,u32> = HashMap::new();
    let mut int_to_minimizer : HashMap<u32,String> = HashMap::new();
    let mut minim_idx : u32 = 0;
    for minim in all_minimizers.iter() {
        minimizer_to_int.insert(minim.to_string(), minim_idx);
        int_to_minimizer.insert(minim_idx,         minim.to_string());
        minim_idx += 1;
    }

    // create a hash table containing (kmers, count)
    // then will keep only those with count > 1
    let mut dbg_nodes : HashMap<Kmer,u32> = HashMap::new(); // it's a Counter
    for read in somewhat_reads
    {
        let read_transformed : Vec<u32> = read.iter().map(|minim| minimizer_to_int[minim]).collect();
        
        if read_transformed.len() <= k { continue; }

        for i in 0..(read_transformed.len()-k+1)
        {
            let mut node : Kmer = Kmer::make_from(&read_transformed[i..i+k]);
            if revcomp_aware { node = node.normalize(); }
            *dbg_nodes.entry(node).or_insert(0) += 1;
        }
    }

    println!("nodes before abund-filter: {}", dbg_nodes.len());
    dbg_nodes.retain(|x,c| c > &mut 1);
    println!("              nodes after: {}", dbg_nodes.len());

    let mut dbg_edges : Vec<(&Kmer,&Kmer)> = Vec::new();

    // index k-1-mers
    let mut km_index : HashMap<Overlap,Vec<&Kmer>> = HashMap::new();
    for node in dbg_nodes.keys()
    {
        let first :Overlap = node.prefix().normalize();
        let second  :Overlap = node.suffix().normalize();
        let mut insert_km = |key,val| 
        {
            match km_index.entry(key) {
                Entry::Vacant(e) => { e.insert(vec![val]); },
                Entry::Occupied(mut e) => { e.get_mut().push(val); }
            }
        };
        insert_km(first,node);
        insert_km(second,node);
    }

    // create a vector of dbg edges
    for n1 in dbg_nodes.keys() 
    {

        /*
        let maybe_insert_edge = |n1 :&[u32;k], n2 :&[u32;k]| {
             if n1[1..] == n2[0..k-1] {
                dbg_edges.push((n1,n2));
             }
        };
        */

        // bit of a rust noob way to code this, because i'm not too familiar with types yet..
        let key1=n1.suffix().normalize();
        let key2=n1.prefix().normalize();
        for key in [key1,key2].iter()
        {
            if km_index.contains_key(&key)
            {
                let list_of_n2s : &Vec<&Kmer> = km_index.get(&key).unwrap();
                for n2 in list_of_n2s {
                    let rev_n2 = n2.reverse();
                    if n1.suffix() == n2.prefix() {
                        dbg_edges.push((n1,n2));
                    }
                    // I wanted to do this closure, but try it, it doesn't work. A borrowed data problem
                    // apparently.
                    //maybe_insert_edge(n1,n2);
                    if revcomp_aware {
                        if n1.prefix() == n2.suffix()
                          || n1.prefix() == rev_n2.suffix()
                          || n1.suffix() == rev_n2.prefix() {
                            dbg_edges.push((n1,n2));
                        }

                    }
                }
            }
        }
    }
    println!("edges: {}", dbg_edges.len());

    // create a real bidirected dbg object using petgraph
    let mut gr = DiGraph::<Kmer,Kmer>::new();
    let mut node_indices : HashMap<Kmer,NodeIndex> = HashMap::new(); // bit redundant info, as nodes indices are in order of elements in dbg_nodes already; but maybe don't want to binary search inside it.
    for node in dbg_nodes.keys() { 
        let index = gr.add_node(node.clone());
        node_indices.insert(node.clone(),index);
    }
    let vec_edges : Vec<(NodeIndex,NodeIndex)> = dbg_edges.iter().map(|(n1,n2)| (node_indices.get(&n1).unwrap().clone(),node_indices.get(&n2).unwrap().clone())).collect();

    gr.extend_with_edges( vec_edges );

    // graphml output
    if false
    //if true
    {
        let graphml = GraphMl::new(&gr).pretty_print(true);
        std::fs::write("graph.graphml", graphml.to_string()).unwrap();
    }
    let output_graph_filename = "graph.gfa";

    // gfa output
    gfa_output::output_gfa(&gr, &dbg_nodes, output_graph_filename);

    // write sequences of minimizers for each node
    // and also read sequences corresponding to those minimizers
    seq_output::write_minimizers_and_seq_of_kmers(output_graph_filename, &node_indices);
}
