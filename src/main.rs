#![allow(unused_variables)]
#![allow(non_upper_case_globals)]

use pbr::ProgressBar;
use bio::io::fasta;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use petgraph_graphml::GraphMl;
use petgraph::graph::DiGraph;
use petgraph::graph::NodeIndex;
//use adler32::RollingAdler32;
use std::fs;
use structopt::StructOpt;
use std::path::PathBuf;
use itertools::iproduct;

mod utils;
mod gfa_output;
mod seq_output;
mod minimizers;

mod kmer_vec;
// mod kmer_array; // not working yet

const revcomp_aware: bool = true; // shouldn't be set to false except for strand-directed data or for debugging

//use typenum::{U31,U32}; // for KmerArray
type Kmer = kmer_vec::KmerVec;
type Overlap= kmer_vec::KmerVec;

pub struct Params
{
    l: usize,
    k: usize,
    percentage_retain_hashes :f64,
    size_miniverse: u32,
    average_lmer_count : f64,
    min_kmer_abundance : usize,
    levenshtein_minimizers : usize,
}

fn extract_minimizers(seq: &str, params: &Params, lmer_counts: &HashMap<String,u32>, minimizer_to_int: &HashMap<String,u32>) -> (Vec<String>, Vec<u32>, Vec<u32>)
{
    minimizers::minhash(seq, params, lmer_counts, minimizer_to_int)
        //wk_minimizers(seq, percentage_retain_hashes) // unfinished
}

fn debug_output_read_minimizers(seq_str: &String, read_minimizers : &Vec<String>, read_minimizers_pos :&Vec<u32>)
{
    println!("seq: {}",seq_str);
    print!("min: ");
    let mut current_minimizer :String = "".to_string();
    for i in 0..seq_str.len()
    {
        if read_minimizers_pos.contains(&(i as u32))
        {
            let index = read_minimizers_pos.iter().position(|&r| r == i as u32).unwrap();
            current_minimizer = read_minimizers[index].clone();
            let c = current_minimizer.remove(0);
            if c == seq_str.chars().nth(i).unwrap()
            {
                print!("X");
            }
                else
            {
                print!("x");
            }
            continue;
        }
        if current_minimizer.len() > 0
        {
            let c = current_minimizer.remove(0);
            print!("{}",c);
        }
        else
        {
            print!(".");
        }
    }
    println!("");

}

// https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter
fn kproduct(seq: String, k: u32) -> Vec<String> {
    match k {
        0 => vec![],
        1 => seq.chars().map(|c| c.to_string()).collect(),
        2 => iproduct!(seq.chars(), seq.chars()).map(|(a, b)| format!("{}{}", a, b)).collect(),
        _ => iproduct!(kproduct(seq.clone(), k - 1), seq.chars()).map(|(a, b)| format!("{}{}", a, b)).collect(),
    }
}



#[derive(Debug, StructOpt)]
#[structopt(name = "rust-mhdbg", about = "Original implementation of MinHash de Bruijn graphs")]
struct Opt {
    /// Activate debug mode
    // short and long flags (-d, --debug) will be deduced from the field's name
    #[structopt(short, long)]
    debug: bool,

    /// Input file
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,

    /// Output graph/sequences prefix 
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,

    #[structopt(short, long)]
    k: Option<usize>,
    #[structopt(short, long)]
    l: Option<usize>,
    #[structopt(long)]
    pch: Option<f64>,
    #[structopt(long)]
    minabund: Option<usize>,
    #[structopt(long)]
    levenshtein_minimizers: Option<usize>,
    #[structopt(long)]
    test1: bool,
    #[structopt(long)]
    test2: bool,
}


fn main() {
    let opt = Opt::from_args();      

    let mut filename = PathBuf::new();
    let mut output_prefix;
    let mut k: usize = 10;
    let mut l: usize = 12;
    let mut percentage_retain_hashes :f64 = 0.10;
    let mut min_kmer_abundance: usize = 2;
    let mut levenshtein_minimizers: usize = 0;

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
    if !opt.minabund.is_none() { min_kmer_abundance = opt.minabund.unwrap() } else { println!("Warning: using default min kmer abundance value ({})",min_kmer_abundance); }

    if !opt.levenshtein_minimizers.is_none() { levenshtein_minimizers = opt.levenshtein_minimizers.unwrap() }
    let minimizer_type = match levenshtein_minimizers { 0 => "reg", 1 => "lev1", 2 => "lev2",_ => "levX" };
    if opt.levenshtein_minimizers.is_none() { println!("Warning: using default minimizer type ({})",minimizer_type); }


    output_prefix = PathBuf::from(format!("{}graph-k{}-p{}-l{}",minimizer_type,k,percentage_retain_hashes,l));

    if !opt.reads.is_none() { filename = opt.reads.unwrap().clone(); } 
    if !opt.prefix.is_none() { output_prefix = opt.prefix.unwrap(); } else { println!("Warning: using default prefix ({})",output_prefix.to_str().unwrap()); }
    
    if filename.as_os_str().is_empty() { panic!("please specify an input file"); }

    let debug = opt.debug;

    let size_miniverse = match revcomp_aware
    {
        false => 4f32.powf(l as f32) as u32,
        true => 4f32.powf(l as f32) as u32 / 2
    };

    let mut params = Params { 
        l,
        k,
        percentage_retain_hashes,
        size_miniverse,
        average_lmer_count: 0.0,
        min_kmer_abundance,
        levenshtein_minimizers,
    };

    // init some useful objects
    let mut nb_minimizers_per_read : f64 = 0.0;
    let mut nb_reads : u64 = 0;
    // get file size for progress bar
    let metadata = fs::metadata(&filename).expect("error opening input file");
    let file_size = metadata.len();
    let mut pb = ProgressBar::new(file_size);

    let (minimizer_to_int, int_to_minimizer, lmer_counts) = minimizers::minimizers_preparation(&mut params, &filename, file_size, levenshtein_minimizers);

    // fasta parsing
    // possibly swap it later for https://github.com/aseyboldt/fastq-rs
    let reader = fasta::Reader::from_file(&filename).unwrap();

    // create a hash table containing (kmers, count)
    // then will keep only those with count > 1
    let mut dbg_nodes   : HashMap<Kmer,u32> = HashMap::new(); // it's a Counter
    let mut kmer_seqs   : HashMap<Kmer,String> = HashMap::new(); // associate a dBG node to its sequence
    let mut minim_shift : HashMap<Kmer,(u32,u32)> = HashMap::new(); // records position of second minimizer in sequence

    for result in reader.records() {
        let record = result.unwrap();
        let seq = record.seq();

        let seq_str = String::from_utf8_lossy(seq);

        //println!("seq: {}", seq_str);

        let (read_minimizers, read_minimizers_pos, read_transformed) = extract_minimizers(&seq_str, &params, &lmer_counts, &minimizer_to_int);

        // stats
        nb_minimizers_per_read += read_minimizers.len() as f64;
        nb_reads += 1;
        pb.add(record.seq().len() as u64 + record.id().len() as u64); // get approx size of entry

        // some debug
        //debug_output_read_minimizers(&seq_str.to_string(), &read_minimizers, &read_minimizers_pos);

        if read_transformed.len() <= k { continue; }

        for i in 0..(read_transformed.len()-k+1)
        {
            let mut node : Kmer = Kmer::make_from(&read_transformed[i..i+k]);
            let mut seq_reversed = false;
            if revcomp_aware { 
                let (node_norm, reversed) = node.normalize(); 
                node = node_norm;
                seq_reversed = reversed;
            }
            let entry = dbg_nodes.entry(node.clone()).or_insert(0);
            *entry += 1;

            // record sequences associated to solid kmers
            if *entry == min_kmer_abundance as u32
            {
                let mut seq = seq_str[read_minimizers_pos[i] as usize..(read_minimizers_pos[i+k-1] as usize + l)].to_string();
                if seq_reversed {
                    seq = utils::revcomp(&seq);
                }
                kmer_seqs.insert(node.clone(), seq.clone());
                let position_of_second_minimizer = match seq_reversed {
                    true => read_minimizers_pos[i+k-1]-read_minimizers_pos[i+k-2],
                    false => read_minimizers_pos[i+1]-read_minimizers_pos[i]
                };
                let position_of_second_to_last_minimizer = match seq_reversed {
                    true => read_minimizers_pos[i+1]-read_minimizers_pos[i],
                    false => read_minimizers_pos[i+k-1]-read_minimizers_pos[i+k-2]
                };

                minim_shift.insert(node.clone(), (position_of_second_minimizer, position_of_second_to_last_minimizer));

                // some sanity checks
                if levenshtein_minimizers == 0
                {
                    for minim in &read_minimizers[i..i+k-1]
                    {
                        debug_assert!((!&seq.find(minim).is_none()) || (!utils::revcomp(&seq).find(minim).is_none()));
                    }
                }
            }
        }        
    }
    pb.finish_print("done converting reads to minimizers");
    nb_minimizers_per_read /= nb_reads as f64;

    println!("avg number of minimizers/read: {}",nb_minimizers_per_read);

    println!("nodes before abund-filter: {}", dbg_nodes.len());
    dbg_nodes.retain(|x,&mut c| c >= (min_kmer_abundance as u32));
    println!("              nodes after: {}", dbg_nodes.len());

    let mut dbg_edges : Vec<(&Kmer,&Kmer)> = Vec::new();

    // index k-1-mers
    let mut km_index : HashMap<Overlap,Vec<&Kmer>> = HashMap::new();
    for node in dbg_nodes.keys()
    {
        let first :Overlap = node.prefix().normalize().0;
        let second  :Overlap = node.suffix().normalize().0;
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
        let rev_n1 = n1.reverse();

        /*
           let maybe_insert_edge = |n1 :&[u32;k], n2 :&[u32;k]| {
           if n1[1..] == n2[0..k-1] {
           dbg_edges.push((n1,n2));
           }
           };
           */

        // bit of a rust noob way to code this, because i'm not too familiar with types yet..
        let key1=n1.suffix().normalize().0;
        let key2=n1.prefix().normalize().0;
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
                        if n1.suffix() == rev_n2.prefix()
                            || rev_n1.suffix() == n2.prefix()
                                || rev_n1.suffix() == rev_n2.prefix() {
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

    // gfa output
    println!("writing GFA..");
    gfa_output::output_gfa(&gr, &dbg_nodes, &output_prefix, &kmer_seqs, &int_to_minimizer, &minim_shift, levenshtein_minimizers);

    // write sequences of minimizers for each node
    // and also read sequences corresponding to those minimizers
    println!("writing sequences..");
    seq_output::write_minimizers_and_seq_of_kmers(&output_prefix, &node_indices, &kmer_seqs, &dbg_nodes, k, l);
}
