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
    let n = params.n;
    let k = params.k;
    let mut bucket_seqs = Vec::<Vec<u32>>::new();
    let mut potential_kmers : HashMap<u32, Vec<Vec<u32>>> = HashMap::new();
    let mut first : bool = true;
    let mut final_cons = Vec::<u32>::new();
    let mut scoring = poa::Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });
    let mut aligner = poa::Aligner::new(scoring, read_transformed.to_vec());
    let mut prev_len = 0;	
    let mut min_prev_len = 99999;
    let mut start_pos = 0;
    let mut offset_hash : HashMap<Vec<u32>, u32> = HashMap::new(); 
    for i in 0..read_transformed.len()-n+1 {
        let mut pileup_seqs = Vec::<Vec<u32>>::new();
        let bucket_idx = read_transformed[i..i+n].to_vec();
        let entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<Vec<u32>>::new());
        //entry.dedup();
        for query in entry.iter() {
            if !bucket_seqs.contains(query) {
                    pileup_seqs.push(query.to_vec());
                    bucket_seqs.push(query.to_vec());
                    let mut offset;	
                    let mut offset_reg = query.iter().position(|&x| x == bucket_idx[0]);	
                    offset = offset_reg.unwrap();	
                    if offset > prev_len {prev_len = offset;}	
                    if offset < min_prev_len {min_prev_len = offset;}    
            }
           // }
        // }
        
        }
        for seq in pileup_seqs.iter() {	
            //print!("Seq\t");	
            //for min in seq.iter() {	
            //    print!("{}\t", min);	
        // }	            // }
            //print!("\n");	
            let mut new_seq = Vec::<u32>::new();	
            let mut offset;	
            let mut offset_reg = seq.iter().position(|&x| x == bucket_idx[0]);	
            offset = offset_reg.unwrap();	
            let mut new_seq = Vec::<u32>::new();	
            if prev_len+i != offset {	
                for _ in offset..prev_len+i {	
                    new_seq.push(0)	
                }	
            }       
            offset_hash.insert(seq.to_vec(), ((prev_len+i) as usize - offset) as u32);	
            for min in seq.iter() {	
                new_seq.push(*min);	
            }     	
           // print!("New\t");	
            /*for min in new_seq.iter() {	
                print!("{}\t", min);	
            }	
            print!("\n");	*/
        }
        pileup_seqs.sort_by_key(|x| offset_hash[&x.to_vec()]);
        for i in 0..pileup_seqs.len() {	
            aligner.global(pileup_seqs[i].to_vec());
            aligner.add_to_graph();
        }
            
        
        
        
        
        //println!("OG\t{:?}", og_kmer);
        //println!("After\t{:?}", consensus);
        //bucket_seqs.push(consensus.to_vec());
        
       // let graphml = GraphMl::new(&aligner.poa.graph).pretty_print(true);
        //std::fs::write("graph.graphml", graphml.to_string()).unwrap();
            //println!("Consensus:\t{}", minimizers::normalize_minimizer(&std::str::from_utf8(&consensus).unwrap().to_string()));
        //}
    }
    //aligner.poa.get_max_path();
    let (mut cons_scores, mut cons_next) = aligner.poa.set_max_weight_scores(&params);
    let mut consensus = aligner.poa.find_consensus(cons_scores, cons_next);
   /* for i in 0..bucket_seqs.len() {
        if i == 0 {
            for lmer in bucket_seqs[i].iter() {
                final_cons.push(*lmer);
            }
        }
        else {
            if !final_cons.contains(&bucket_seqs[i][bucket_seqs[i].len()-1]) {
                final_cons.push(bucket_seqs[i][bucket_seqs[i].len()-1]);
            }
        }
    }*/
    println!("OG\t{:?}", read_transformed);
    println!("After\t{:?}", consensus);
   // let (mut cons_scores, mut cons_next) = aligner.poa.set_max_weight_scores();
    //let mut consensus = aligner.poa.find_consensus(cons_scores, cons_next);
   // println!("OG\t{:?}", read_transformed);
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
