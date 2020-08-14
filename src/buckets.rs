use std::collections::HashMap;
use super::Kmer;
use array_tool::vec::Intersect;
use super::kmer_vec;
use super::Params;
use super::utils;
use std::collections::HashSet;
use super::minimizers;
use super::extract_minimizers;
use rust_spoa::poa_consensus;
use petgraph_graphml::GraphMl;
use petgraph::graph::EdgeIndex;
use std::io::BufWriter;
use super::poa;
use super::poa_new;
use super::jaccard_distance;
use super::ec_reads;
use std::fs::File;


pub fn query_buckets(pairwise_jaccard : &mut HashMap<(Vec<u64>, Vec<u64>), f64>, ec_file_poa: &mut BufWriter<File>, read_ids : &mut HashMap<Vec<u64>, String>, mut corrected : &mut HashMap<Vec<u64>, Vec<u64>>, read_transformed : Vec<u64>, buckets : &mut HashMap<Vec<u64>, Vec<Vec<u64>>>, params : &Params) -> Vec<u64>{
    let n = params.n;
    let k = params.k;
    let mut bucket_seqs = Vec::<Vec<u64>>::new();
    let mut scoring = poa::Scoring::new(-1, 0, |a: u64, b: u64| if a == b { 1i32 } else { -1i32 });
    let mut aligner = poa::Aligner::new(scoring, read_transformed.to_vec());
    let mut aligned : HashMap<Vec<u64>, bool> = HashMap::new();
    let mut poa_ids = Vec::<String>::new();
    let mut seq_id = read_ids[&read_transformed.to_vec()].to_string();
    for i in 0..read_transformed.len()-n+1 {
        let bucket_idx = read_transformed[i..i+n].to_vec();
        let entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<Vec<u64>>::new());
        //entry.dedup();
        for query in entry.iter() {
            if query.to_vec() == read_transformed.to_vec() {continue;}
            let entry = aligned.entry(query.to_vec()).or_insert(false);
            if !*entry {
                poa_ids.push(read_ids[&query.to_vec()].to_string());
                let tuple = (read_transformed.to_vec(), query.to_vec());
                let tuple_rev = (tuple.1.to_vec(), tuple.0.to_vec());
                let mut similarity = 0.0;
                if pairwise_jaccard.contains_key(&tuple) {
                    similarity = pairwise_jaccard[&tuple];
                }
                else if pairwise_jaccard.contains_key(&tuple_rev) {
                    similarity = pairwise_jaccard[&tuple_rev];
                }
                else {
                    similarity = jaccard_distance(read_transformed.to_vec(), query.to_vec());
                    pairwise_jaccard.insert(tuple, similarity);
                    pairwise_jaccard.insert(tuple_rev, similarity);
                }
                *entry = true;   
                if similarity < 0.33 {
                    //println!("Jaccard Similarity: {}", similarity);
                    continue;
                }

                //if common.len() <= n+1 {
                //    continue;
                //}
                    //let aln = poa_graph.align_sequence(query.to_vec());
                    //poa_graph.incorporate_alignment(aln, &vec![1], query.to_vec());
                    //consensus = poa_graph.consensus();
                aligner.global(query.to_vec());
                aligner.add_to_graph();
                if similarity >= 0.33 {
                    bucket_seqs.push(query.to_vec());
                }
                //let mut offset;	
                //let mut offset_reg = query.iter().position(|&x| x == bucket_idx[0]);	
                //offset = offset_reg.unwrap();	
                //if offset > prev_len {prev_len = offset;}	
                //if offset < min_prev_len {min_prev_len = offset;} 
            }
           // }
        // }
        
        }
    }
    let mut consensus = aligner.poa.consensus();
    for seq in bucket_seqs.iter() {
        *corrected.entry(seq.to_vec()).or_insert(Vec::<u64>::new()) = consensus.to_vec();
    }
    ec_reads::record_poa(ec_file_poa, &seq_id, poa_ids);
    consensus.to_vec()
    
}
    
pub fn buckets_insert(seq : Vec<u64>, i : usize, buckets : &mut HashMap<Vec<u64>, Vec<Vec<u64>>>, mut dbg_nodes : &mut HashMap<Kmer, u32>) { 
    for j in 0..seq.len()-i+1 {
        let mut bucket_idx = Vec::<u64>::new();
        bucket_idx = seq[j..i+j].to_vec();
        let mut entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<Vec<u64>>::new());
       // for _ in 0..sub_counts[&bucket_idx] {
            entry.push(seq.to_vec());
        //}

    }
    

}
