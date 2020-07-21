use std::collections::HashMap;
use super::Kmer;
use super::kmer_vec;
use super::Params;
use super::utils;
use super::extract_minimizers;

pub fn get_kmers(dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> Vec<Vec<u32>> {
    dbg_nodes.iter().map(|(kmer, count)| kmer_vec::get(kmer).to_vec()).collect()
}

pub fn enumerate_buckets(mut dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>> {
    let kmers = get_kmers(&mut dbg_nodes, &params);
    let mut buckets : HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>> = HashMap::new();
    for kmer in kmers.iter() {
        buckets_insert(kmer.to_vec(), 2, &mut buckets, &mut dbg_nodes);
    }
    for (key, entry) in &mut buckets {
        entry.sort_by_key(|tuple| tuple.0.iter().position(|&x| x == key[0]).unwrap());
        entry.reverse();

    }
    println!("{:?}", buckets.len());
    //buckets.retain(|key, entry| entry.len() > 3 && entry.iter().filter(|tuple| tuple.1 == 1).count() == entry.len());
    println!("{:?}", buckets.len());
    buckets
}
pub fn recur_multiple (consensus: Vec<u32>,  entry : Vec<(Vec<u32>, u32)>, mut multiple : &mut Vec::<Vec<u32>>, counts : &HashMap<u32, u32>) -> Vec::<Vec<u32>>{
    let consensus_suffix = &consensus[1..consensus.len()-1];
    multiple.push(consensus.to_vec());
    for (kmer, count) in entry.iter() {
        if &kmer[0..kmer.len()-1] == consensus_suffix && counts[&kmer[kmer.len()-1]] > 1 {
            recur_multiple(kmer.to_vec(), entry.to_vec(), &mut multiple, &counts);

        }
    }
    multiple.to_vec()


}
pub fn query_buckets(read_transformed : Vec<u32>, mut dbg_nodes: &mut HashMap<Kmer,u32>, buckets : &mut HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>>, mut seq_str: &mut String, params : &Params, mut kmer_seqs_tot: &mut HashMap<Kmer,String>, lmer_counts: &HashMap<String,u32>, minimizer_to_int : &HashMap<String,u32>, int_to_minimizer : &HashMap<u32,String>) -> String{
    let sub = 2;
    let (read_minimizers, read_minimizers_pos, read_transformed) = extract_minimizers(&seq_str, &params, &lmer_counts, &minimizer_to_int);
    let k = params.k;
    let mut consensus_strs = Vec::<String>::new();
    for i in 0..read_transformed.len()-k+1 {
        let og_kmer = &read_transformed[i..i+k];
        let kmer_node = Kmer::make_from(&og_kmer.to_vec());
        let (node_norm, reversed) = kmer_node.normalize();
        let mut counter = 1;
        //println!("OG Kmer {:?}", og_kmer);
        //println!("Seq_str {}", seq_str);
        //println!("Kmer seqs tot {}", kmer_seqs_tot[&node_norm]);
        for j in 0..og_kmer.len()-sub+1 {
            let mut bucket_idx = Vec::<u32>::new();
            for k in 0..sub {
                bucket_idx.push(og_kmer.to_vec()[j+k]);
            }
            let mut entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<(Vec<u32>, u32)>::new());
            if entry.len() == 0 {continue;}
            
            let mut base_counts : HashMap<u32, Vec<String>> = HashMap::new();
            let mut bucket_seqs : HashMap<Vec<u32>, String> = HashMap::new();
            let mut bucket_gapped_seqs : HashMap<Vec<u32>, String> = HashMap::new();
            let mut counts : HashMap<u32, u32> = HashMap::new();
            let mut found : bool = false;
            for (kmer, count) in entry.iter() {
                if kmer.to_vec() == og_kmer.to_vec() {found = true;}
                let kmer_node = Kmer::make_from(&kmer.to_vec());
                let (node_norm, reversed) = kmer_node.normalize();
                //println!("{:?}\t{}\t{:?}\t{:?}", key, dbg_nodes[&kmer_node], kmer, kmer_seqs_tot[&kmer_node]);
                //if reversed {
                 //   bucket_seqs.insert(kmer.to_vec(), utils::revcomp(&kmer_seqs_tot[&node_norm]).to_string());
                //}
               // else {
                    bucket_seqs.insert(kmer.to_vec(), kmer_seqs_tot[&node_norm].to_string());
                //}
            }
            if !found {break;}
           // println!("Found false kmer");
            let mut prev_len = 0;
            let mut min_prev_len = 99999;
            let mut og_str_idx = Vec::<u32>::new();
            for (kmer, seq) in &bucket_seqs {
                if seq == seq_str {og_str_idx = kmer.to_vec()};
                //println!("{}", seq.len());
                let mut offset;
                let mut offset_reg = seq.find(&int_to_minimizer[&bucket_idx[0]]);
                if offset_reg == None {offset = seq.find(&utils::revcomp(&int_to_minimizer[&bucket_idx[0]])).unwrap();}
                else {offset = offset_reg.unwrap();}
                if offset > prev_len {prev_len = offset;}
                if offset < min_prev_len {min_prev_len = offset;}       
            }
            for (kmer, seq) in &bucket_seqs {
                    let mut offset;
                    let mut offset_reg = seq.find(&int_to_minimizer[&bucket_idx[0]]);
                    if offset_reg == None {offset = seq.find(&utils::revcomp(&int_to_minimizer[&bucket_idx[0]])).unwrap();}
                    else {offset = offset_reg.unwrap();}
                let mut new_str = String::new();
                if prev_len != offset {
                    for i in offset..prev_len {
                        new_str.push_str("-")
                    }
                }            
                new_str.push_str(&seq.to_string());
                //println!("{}", new_str);
                bucket_gapped_seqs.insert(kmer.to_vec(), new_str);
            }
            for (kmer, new_str) in &bucket_gapped_seqs {
                for j in 0..new_str.len() {
                    let mut pos_entry = base_counts.entry(j as u32).or_insert(Vec::<String>::new());
                    pos_entry.push(new_str.chars().nth(j).unwrap().to_string());
                }     
            }
            let mut consensus_str = String::new();
            for pos in 0..base_counts.keys().len() {
                let mut base_freq : HashMap<String, u32> = HashMap::new();
                for base in base_counts[&(pos as u32)].iter() {
                    if base == "-" {continue;}
                    let base_cnt = base_freq.entry(base.to_string()).or_insert(0);
                    *base_cnt += 1;
                }
                let max_base = base_freq.into_iter().max_by_key(|&(_, count)| count).map(|(val, _)| val).unwrap();
                consensus_str.push_str(&max_base.clone());
            }
            //println!("For this kmer, consensus {}: {}", counter, consensus_str);
            let mut seq_str_new = String::new();
            for (kmer, new_str) in &bucket_gapped_seqs {
                if kmer.to_vec() == og_str_idx.to_vec() {
                    //println!("here");
                    let mut start_pos = 0;
                    for i in 0..new_str.len() {
                        if new_str.chars().nth(i).unwrap().to_string() != "-" {
                            start_pos = i;
                            break;
                        }
                    }
                    seq_str_new = consensus_str[start_pos..new_str.len()].to_string();
                    break;  
                }
            }
            //println!("Seq_str was");
            //println!("{}", seq_str);
            //println!("Now");
            //println!("{}", seq_str_new.to_string());
            consensus_strs.push(seq_str_new);
            counter += 1;
        }
        break;
    }
    let mut seq_base_counts : HashMap<u32, Vec<String>> = HashMap::new();
    let mut consensus_seq = String::new();
    //println!("Seq str \t {}", seq_str);
    for cand_seq in consensus_strs.iter() {
        //println!("Candidate seq \t {}", cand_seq);
        for j in 0..cand_seq.len() {
            let mut pos_entry = seq_base_counts.entry(j as u32).or_insert(Vec::<String>::new());
            pos_entry.push(cand_seq.chars().nth(j).unwrap().to_string());
        }     
    }
    let mut consensus_seq = String::new();
    for pos in 0..seq_base_counts.keys().len() {
        let mut base_freq : HashMap<String, u32> = HashMap::new();
        for base in seq_base_counts[&(pos as u32)].iter() {
            let base_cnt = base_freq.entry(base.to_string()).or_insert(0);
            *base_cnt += 1;
        }
        let max_base = base_freq.into_iter().max_by_key(|&(_, count)| count).map(|(val, _)| val).unwrap();
        consensus_seq.push_str(&max_base.clone());
    }
    if consensus_seq != "" {
        //println!("Consensus seq \t {}", consensus_seq);
        consensus_seq.to_string()
    }
    else {
        seq_str.to_string()
    }

}

pub fn buckets_insert(kmer : Vec<u32>, i : usize, buckets : &mut HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>>, mut dbg_nodes: &mut HashMap<Kmer,u32>) {
    for j in 0..kmer.len()-i+1 {
        let mut bucket_idx = Vec::<u32>::new();
        for k in 0..i {
            bucket_idx.push(kmer.to_vec()[j+k]);
        }
        let mut entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<(Vec<u32>, u32)>::new());
        let kmer_node = Kmer::make_from(&kmer.to_vec());
        let (node_norm, reversed) = kmer_node.normalize();
        entry.push((kmer.to_vec(), dbg_nodes[&node_norm]));
    }

}