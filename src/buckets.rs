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



pub fn query_buckets(int_to_minimizer: &mut HashMap<u64, String>, read_to_seq_pos: &mut HashMap<Vec<u64>, (String, Vec<u32>)>, pairwise_jaccard : &mut HashMap<(Vec<u64>, Vec<u64>), (f64, f64)>, ec_file_poa: &mut BufWriter<File>, read_ids : &mut HashMap<Vec<u64>, String>, mut corrected : &mut HashMap<Vec<u64>, (String, Vec<String>, Vec<u32>, Vec<u64>)>, read_transformed : &Vec<u64>, buckets : &mut HashMap<Vec<u64>, Vec<Vec<u64>>>, params : &Params) -> (String, Vec<String>, Vec<u32>, Vec<u64>) {
    let n = params.n;
    let k = params.k;
    let l = params.l;
    let mut bucket_seqs = Vec::<&Vec<u64>>::new();
    let mut scoring = poa::Scoring::new(-1, 0, |a: u64, b: u64| if a == b { 1i32 } else { -1i32 }).xclip(0).yclip(0);
    let mut aligner = poa::Aligner::new(scoring, read_transformed);
    let mut aligned : HashMap<&Vec<u64>, bool> = HashMap::new();
    let mut poa_ids = Vec::<String>::new();
    let mut seq_id = read_ids[read_transformed].to_string();
    let mut pair_map : HashMap<(u64, u64), String> = HashMap::new();
    let seq_str = read_to_seq_pos[&read_transformed.to_vec()].0.to_string();
    let read_minimizers_pos = read_to_seq_pos[&read_transformed.to_vec()].1.to_vec();
    for i in 0..read_transformed.len()-1 {
        pair_map.insert((read_transformed[i], read_transformed[i+1]), seq_str[read_minimizers_pos[i] as usize..read_minimizers_pos[i+1] as usize].to_string());
        //println!("{:?}", ((read_transformed[i], read_transformed[i+1]), seq_str[read_minimizers_pos[i] as usize ..read_minimizers_pos[i+1] as usize].to_string()));

    }
    for i in 0..read_transformed.len()-n+1 {
        let mut max_seqs = 0;
        let bucket_idx = read_transformed[i..i+n].to_vec();
        let entry = &buckets[&bucket_idx];
        //entry.dedup();
        for query in entry.iter() {
            if query == read_transformed {continue;}
            let entry = aligned.entry(query).or_insert(false);
            if !*entry {
                bucket_seqs.push(query);
                *entry = true;  



                //if common.len() <= n+1 {
                //    continue;
                //}
                    //let aln = poa_graph.align_sequence(query.to_vec());
                    //poa_graph.incorporate_alignment(aln, &vec![1], query.to_vec());
                    //consensus = poa_graph.consensus();
                
                //if containment > 0.6 {
                //    bucket_seqs.push(query);
                //}
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
    let mut new = bucket_seqs.iter().map(|&seq| (seq.to_vec(), jaccard_distance(seq, &read_transformed))).filter(|tuple| ((tuple.1).0 > 0.15 && (tuple.1).1 > 0.25)).collect::<Vec<(Vec<u64>, (f64, f64))>>();
    new.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    new.reverse();
    for seq in new.iter() {
        let query = seq.0.to_vec();
        let seq_str = read_to_seq_pos[&query.to_vec()].0.to_string();
        let read_minimizers_pos = read_to_seq_pos[&query.to_vec()].1.to_vec();
        for i in 0..query.len()-1 {
            //println!("{:?}", ((query[i], query[i+1]), seq_str[read_minimizers_pos[i] as usize ..read_minimizers_pos[i+1] as usize].to_string()));
            pair_map.insert((query[i], query[i+1]), seq_str[read_minimizers_pos[i] as usize ..read_minimizers_pos[i+1] as usize].to_string());
        }
    }
    let mut max_len = 200;
    if new.len() < max_len {max_len = new.len();}
    for i in 0..max_len {
        poa_ids.push(read_ids[&new[i].0].to_string());
       //println!("Containment: {}", (new[i].1).0);
        //println!("Similarity: {}", (new[i].1).1);
        aligner.global(&new[i].0);
        aligner.add_to_graph();
    }
    

    let mut consensus = aligner.poa.consensus(&params);
    
    let consensus_read = consensus.iter().map(|minim| int_to_minimizer[minim].to_string()).collect::<Vec<String>>();
    let mut consensus_str = String::new();
    let mut pos_idx = 0;
    let mut consensus_pos = Vec::<u32>::new();
    if consensus.len() == 0 {return (consensus_str, consensus_read, consensus_pos, consensus)}
    //println!("{}", pair_map.len());
    for i in 0..consensus.len()-1 {
        consensus_pos.push(pos_idx as u32);
        let mut insert = String::new();
        //println!("\n{:?}", (consensus[i], consensus[i+1]));
        if !pair_map.contains_key(&(consensus[i], consensus[i+1])) {
            let alts : Vec<String> = pair_map.clone().into_iter().filter(|(a, b)| a.0 == consensus[i]).map(|tuple| tuple.1).collect();
            //println!("{}", alts.len());
            insert = int_to_minimizer[&consensus[i]].to_string();
            for _ in 0..alts[0].len() {insert.push_str("N")};
        }
        else { insert = pair_map[&(consensus[i], consensus[i+1])].to_string();}
        consensus_str.push_str(&insert);
        pos_idx += insert.len();
    }
    consensus_pos.push(pos_idx as u32);
    consensus_str.push_str(&int_to_minimizer[&consensus[consensus.len()-1]]);
    for (vec, tuple) in new.iter() {
        *corrected.entry(vec.to_vec()).or_insert((String::new(), Vec::<String>::new(), Vec::<u32>::new(), Vec::<u64>::new())) = (consensus_str.to_string(), consensus_read.to_vec(), consensus_pos.to_vec(), consensus.to_vec());
    }
    ec_reads::record_poa(ec_file_poa, &seq_id, poa_ids);
    //println!("{} {} {} {}", consensus_str.len(), consensus_read.len(), consensus_pos.len(), consensus.len());
    (consensus_str, consensus_read, consensus_pos, consensus)
}
    
pub fn buckets_insert(seq : Vec<u64>, i : usize, buckets : &mut HashMap<Vec<u64>, Vec<Vec<u64>>>, mut dbg_nodes : &mut HashMap<Kmer, u32>) { 
    for j in 0..seq.len()-i+1 {
        let mut entry = buckets.entry(seq[j..i+j].to_vec()).or_insert(Vec::<Vec<u64>>::new());
       // for _ in 0..sub_counts[&bucket_idx] {
            entry.push(seq.to_vec());
        //}

    }
    

}

