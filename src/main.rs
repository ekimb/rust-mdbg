#![allow(unused_variables)]
#![allow(non_upper_case_globals)]

use pbr::ProgressBar;
use bio::io::fasta;
use std::collections::HashSet;
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use fasthash::city;
use std::convert::TryInto;
use petgraph_graphml::GraphMl;
use petgraph::graph::DiGraph;
use petgraph::graph::NodeIndex;
use adler32::RollingAdler32;
use std::fs;
use arrayvec::ArrayVec;

mod utils;
mod gfa_output;

const k :usize = 32;
const l :usize = 12;
const percentage_retain_hashes :f32 = 0.01;
const revcomp_aware: bool = true;
 
fn extract_minimizers(seq: &str) -> Vec<String>
{
    minhash(seq)
    //wk_minimizers(seq) // unfinished
}

fn minhash(seq: &str) -> Vec<String>
{
    let hash_mode = 0; // 1: rolling (currently incompat with revcomp), 0: cityhash
    let mut res = Vec::new();
    if seq.len() < l { return res; }

    let mut nl =  4f32.powf(l as f32) as u32;
    if revcomp_aware { nl = nl / 2; }
    let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
    let mut h0 :u32 = 0;
    for start in 0..(&seq.len()-l+1)
    {
        let mut lmer = String::from(&seq[start..start+l]);
        if revcomp_aware {
            let lmer_rev = utils::revcomp(&lmer);
            lmer = std::cmp::min(lmer, lmer_rev);
        }

        if hash_mode == 0 {
            h0 = city::hash32(&lmer) % nl;
        }
        else if hash_mode == 1 {
            if start > 0
            {
                h1.remove(l, seq.as_bytes()[start-1]);
                h1.update(seq.as_bytes()[start + l-1]);
            }
        }
       
        let considered_hash = if hash_mode == 0 { h0 }  else { h1.hash() % nl };
        if considered_hash < (nl as f32*percentage_retain_hashes) as u32 {
            &res.push(lmer.to_string());
            //println!("selected lmer: {} (hash {})", lmer.to_string(), considered_hash);
        }
    }
    res
}

// TODO finish implementing it
#[allow(dead_code)]
fn wk_minimizers(seq: &str) -> Vec<String>
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
    }

    // TODO convert the list of (minimizers,positions) to a vector of list minimizers
    res
}

pub fn reverse_node(node :[u32;k]) -> [u32;k]
{
    let rev_node_tmp :ArrayVec<[u32;k]> = node.into_iter().rev().map(|x| *x).collect();
    rev_node_tmp.into_inner().unwrap()
}

pub fn reverse_overlap(node :[u32;k-1]) -> [u32;k-1]
{
    let rev_node_tmp :ArrayVec<[u32;k-1]> = node.into_iter().rev().map(|x| *x).collect();
    rev_node_tmp.into_inner().unwrap()
}


fn normalize_node(node :[u32;k]) -> [u32;k]
{
    if !revcomp_aware { return node }
    std::cmp::min(node,reverse_node(node))
}

// TODO factorize with previous function
fn normalize_overlap(node :[u32;k-1]) -> [u32;k-1]
{
    if !revcomp_aware { return node }
    std::cmp::min(node,reverse_overlap(node))
}



fn main() {
       
    let mut somewhat_reads : Vec<Vec<String>> = vec![];
    let mut all_minimizers : HashSet<String> = HashSet::new();
    let mut nb_minimizers_per_read : f64 = 0.0;
    let mut nb_reads : u64 = 0;

    let filename = "../read50x_ref10K_e001.fa"; // todo argument
    let filename = "../SRR9969842_vs_chr4.fasta"; // todo argument

    // get file size for progress bar
    let metadata = fs::metadata(filename).expect("error opening input file");
    let file_size = metadata.len();
    let mut pb = ProgressBar::new(file_size);

    // fasta parsing
    // possibly swap it later for https://github.com/aseyboldt/fastq-rs
    let reader = fasta::Reader::from_file(filename).unwrap();
    for result in reader.records() {
	    let record = result.unwrap();
	    let seq = record.seq();
	
        let seq_str = String::from_utf8_lossy(seq);

        //println!("seq: {}", seq_str);

        let minimizers = extract_minimizers(&seq_str);

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

    // create dbg nodes
    let mut dbg_nodes : HashMap<[u32;k],u32> = HashMap::new(); // it's a Counter
    for read in somewhat_reads
    {
        let read_transformed : Vec<u32> = read.iter().map(|minim| minimizer_to_int[minim]).collect();
        
        if read_transformed.len() <= k { continue; }

        for i in 0..(read_transformed.len()-k+1)
        {
            let mut node : [u32;k] = read_transformed[i..i+k].try_into().unwrap();
            if revcomp_aware { node = normalize_node(node); }
            *dbg_nodes.entry(node).or_insert(0) += 1;
        }
    }

    println!("nodes before abund-filter: {}", dbg_nodes.len());
    dbg_nodes.retain(|&x,c| c > &mut 1);
    println!("              nodes after: {}", dbg_nodes.len());

    let mut dbg_edges : Vec<(&[u32;k],&[u32;k])> = Vec::new();

    // index k-1-mers
    let mut km_index : HashMap<[u32;k-1],Vec<[u32;k]>> = HashMap::new();

    for &node in dbg_nodes.keys()
    {
        let first :[u32;k-1]= normalize_overlap(node[0..(k-1)].try_into().unwrap());
        let second  :[u32;k-1]= normalize_overlap(node[1..k].try_into().unwrap());
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
        let key1=normalize_overlap(n1[1..k].try_into().unwrap());
        let key2=normalize_overlap(n1[0..k-1].try_into().unwrap());
        for &key in [key1,key2].iter()
        {
            if km_index.contains_key(&key)
            {
                let list_of_n2s : &Vec<[u32;k]> = km_index.get(&key).unwrap();
                for n2 in list_of_n2s {
                    let rev_n2_tmp :ArrayVec<[u32;k]> = n2.into_iter().rev().map(|x| *x).collect();
                    let rev_n2 = rev_n2_tmp.into_inner().unwrap();
                    if n1[1..] == n2[0..k-1] {
                        dbg_edges.push((n1,n2));
                    }
                    // I wanted to do this closure, but try it, it doesn't work. A borrowed data problem
                    // apparently.
                    //maybe_insert_edge(n1,n2);
                    if revcomp_aware {
                        if n1[..k-1] == n2[1..k] 
                          || n1[..k-1] == rev_n2[1..k]  
                          || n1[1..k] == rev_n2[..k-1] {
                            dbg_edges.push((n1,n2));
                        }

                    }
                }
            }
        }
    }
    println!("edges: {}", dbg_edges.len());


    // create a real graph object
    let mut gr = DiGraph::<[u32;k],[u32;k]>::new();
    let mut node_indices : HashMap<[u32;k],NodeIndex> = HashMap::new(); 
    for &node in dbg_nodes.keys() { 
        let index = gr.add_node(node);
        node_indices.insert(node,index);
    }
    let vec_edges : Vec<(NodeIndex,NodeIndex)> = dbg_edges.iter().map(|(&n1,&n2)| (node_indices.get(&n1).unwrap().clone(),node_indices.get(&n2).unwrap().clone())).collect();

    gr.extend_with_edges( vec_edges );

    // graphml output
    if false
    //if true
    {
        let graphml = GraphMl::new(&gr).pretty_print(true);
        std::fs::write("graph.graphml", graphml.to_string()).unwrap();
    }

    // gfa output
    gfa_output::output_gfa(&gr, &dbg_nodes);
}
