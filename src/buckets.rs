use std::collections::HashMap;
use super::Kmer;
use super::kmer_vec;
use super::Params;
use super::utils;
use super::minimizers;
use super::extract_minimizers;
use rust_spoa::poa_consensus;
use petgraph_graphml::GraphMl;
use super::poa;


fn is_sub<T: PartialEq>(haystack: &[T], needle: &[T]) -> bool {
    haystack.windows(needle.len()).any(|c| c == needle)
}

pub fn get_kmers(dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> Vec<Vec<u32>> {
    dbg_nodes.iter().map(|(kmer, count)| kmer_vec::get(kmer).to_vec()).collect()
}

pub fn enumerate_buckets(mut seq_mins : &mut Vec<Vec<u32>>, mut dbg_nodes: &mut HashMap<Kmer,u32>, mut sub_counts: &mut HashMap<Vec<u32>,u32>, kmer_seqs_tot : &HashMap<Kmer, String>, params : &Params) -> HashMap<Vec<u32>, Vec<Vec<u32>>> {
    let mut buckets : HashMap<Vec<u32>, Vec<Vec<u32>>> = HashMap::new();
    for seq in seq_mins.iter() {
            buckets_insert(seq.to_vec(), params.n, &mut buckets, &mut dbg_nodes, &mut sub_counts, kmer_seqs_tot);
    }
    //for (key, entry) in &mut buckets {
      //  entry.sort_by_key(|seq| seq.iter().position(|&x| x == key[0]).unwrap());
        //entry.reverse();

    //}

    println!("{:?}", buckets.len());
    //buckets.retain(|key, entry| entry.len() > 2);
    println!("{:?}", buckets.len());
    buckets
}

pub fn query_buckets(mut kmer_seqs_tot : &mut HashMap<Kmer,String>, read_transformed : Vec<u32>, mut sub_counts: &mut HashMap<Vec<u32>,u32>, buckets : &mut HashMap<Vec<u32>, Vec<Vec<u32>>>, mut seq_str: &mut String, params : &Params, mut seq_mins: &mut Vec<Vec<u32>>, lmer_counts: &HashMap<String,u32>, minimizer_to_int : &HashMap<String,u32>, int_to_minimizer : &HashMap<u32,String>) -> Vec<u32>{
    let sub = params.n;
    let k = params.k;
    let mut bucket_seqs = Vec::<Vec<u32>>::new();
    let mut potential_kmers : HashMap<u32, Vec<Vec<u32>>> = HashMap::new();
    let mut first : bool = true;
    let mut scoring = poa::Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });
    let mut aligner = poa::Aligner::new(scoring, read_transformed.to_vec());;
    for i in 0..read_transformed.len()-sub+1 {
        let bucket_idx = read_transformed[i..i+sub].to_vec();
        let entry = buckets.entry(bucket_idx).or_insert(Vec::<Vec<u32>>::new());
        for seq in entry.iter() {
                //if !bucket_seqs.contains(seq) {
                    aligner.global(seq.to_vec());
                    aligner.add_to_graph();
                    bucket_seqs.push(seq.to_vec());
               // }
            // }
            
        }
        
       // let graphml = GraphMl::new(&aligner.poa.graph).pretty_print(true);
        //std::fs::write("graph.graphml", graphml.to_string()).unwrap();
            //println!("Consensus:\t{}", minimizers::normalize_minimizer(&std::str::from_utf8(&consensus).unwrap().to_string()));
        //}
    }
    let (mut cons_scores, mut cons_next) = aligner.poa.set_max_weight_scores();
    let mut consensus = aligner.poa.find_consensus(cons_scores, cons_next);
    //println!("OG\t{:?}", read_transformed);
    //println!("After\t{:?}", consensus);
    /*let mut start_pos = 0;
    let mut end_pos = consensus.len();
    
    for i in 0..read_transformed.len() {
        let pos = consensus.iter().position(|element| element == &read_transformed[read_transformed.len()-i-1]);
        if pos == None {continue;}
        else {
            end_pos = pos.unwrap()+1;
            break;
        }
    }
    if end_pos >= read_transformed.len() {start_pos = end_pos - read_transformed.len();}
    let new_consensus = &consensus[start_pos..end_pos];
    println!("Before\t{:?}", consensus);
    println!("After\t{:?}", new_consensus);

    new_consensus.to_vec()*/
    consensus
    
}
pub fn buckets_insert(seq : Vec<u32>, i : usize, buckets : &mut HashMap<Vec<u32>, Vec<Vec<u32>>>, mut dbg_nodes : &mut HashMap<Kmer, u32>, mut sub_counts : &mut HashMap<Vec<u32>, u32>, kmer_seqs_tot : &HashMap<Kmer, String>) {
    /*let mut node : Kmer = Kmer::make_from(&seq);
    let mut seq_reversed = false;
    let (node_norm, reversed) = node.normalize(); 
    node = node_norm;
    seq_reversed = reversed;
    let entry_num = dbg_nodes_all.entry(node.clone()).or_insert(0);
    *entry_num = dbg_nodes[&node];*/
    
    for j in 0..seq.len()-i+1 {
        let mut bucket_idx = Vec::<u32>::new();
        bucket_idx = seq[j..i+j].to_vec();
        let mut entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<Vec<u32>>::new());
       // for _ in 0..sub_counts[&bucket_idx] {
            entry.push(seq.to_vec());
        //}

    }
    

}
