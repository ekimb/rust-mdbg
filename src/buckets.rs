use std::collections::HashMap;
use super::Kmer;
use super::kmer_vec;
use super::Params;
use super::utils;

pub fn get_kmers(dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> Vec<Vec<u32>> {
    dbg_nodes.iter().map(|(kmer, count)| kmer_vec::get(kmer).to_vec()).collect()
}

pub fn enumerate_buckets(mut dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params, kmer_seqs_tot: &mut HashMap<Kmer,String>, minim_shift : &mut HashMap<Kmer,(u32,u32)>) -> HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>> {
    let kmers = get_kmers(&mut dbg_nodes, &params);
    let mut buckets : HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>> = HashMap::new();
    for kmer in kmers.iter() {
        buckets_insert(kmer.to_vec(), 9, &mut buckets, &mut dbg_nodes);
    }
    for (key, entry) in &mut buckets {
        entry.sort_by_key(|tuple| tuple.0.iter().position(|&x| x == key[0]).unwrap());
        entry.reverse();

    }
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
pub fn get_consensus(mut dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params, kmer_seqs_tot: &mut HashMap<Kmer,String>, minim_shift : &mut HashMap<Kmer,(u32,u32)>, int_to_minimizer : &HashMap<u32,String>) -> HashMap<Vec<u32>, String> {
    let mut buckets = enumerate_buckets(&mut dbg_nodes, &params, kmer_seqs_tot, minim_shift);
    println!("{:?}", buckets.len());
    buckets.retain(|key, entry| entry.len() >= 100 && entry.iter().filter(|tuple| tuple.1 == 1).count() == entry.len());
    println!("{:?}", buckets.len());
    let mut seq_switch : HashMap<Vec<u32>, String> = HashMap::new();
    let mut counter = 0; 
    for (key, entry) in &buckets {
        println!("Entry {}", counter);
        let mut base_counts : HashMap<u32, Vec<String>> = HashMap::new();
        let mut bucket_seqs : HashMap<Vec<u32>, String> = HashMap::new();
        let mut bucket_gapped_seqs : HashMap<Vec<u32>, String> = HashMap::new();
        let mut max_index = -1;
        let mut min_index = params.k+1;
        let mut entry_collect = Vec::<u32>::new();
        let mut counts : HashMap<u32, u32> = HashMap::new();
        for (kmer, count) in entry.iter() {
            let kmer_node = Kmer::make_from(&kmer.to_vec());
            let (node_norm, reversed) = kmer_node.normalize();
            //println!("{:?}\t{}\t{:?}\t{:?}", key, dbg_nodes[&kmer_node], kmer, kmer_seqs_tot[&kmer_node]);
            if reversed {
                bucket_seqs.insert(kmer.to_vec(), utils::revcomp(&kmer_seqs_tot[&kmer_node]).to_string());
            }
            else {
                bucket_seqs.insert(kmer.to_vec(), kmer_seqs_tot[&kmer_node].to_string());
            }
        }
        let mut prev_len = 0;
        let mut min_prev_len = 99999;
        for (kmer, seq) in &bucket_seqs {
            //println!("{}", seq.len());
            let mut offset;
            let mut offset_reg = seq.find(&int_to_minimizer[&key[0]]);
            if offset_reg == None {offset = seq.find(&utils::revcomp(&int_to_minimizer[&key[0]])).unwrap();}
            else {offset = offset_reg.unwrap();}
            if offset > prev_len {prev_len = offset;}
            if offset < min_prev_len {min_prev_len = offset;}       
        }
        for (kmer, seq) in &bucket_seqs {
                let mut offset;
                let mut offset_reg = seq.find(&int_to_minimizer[&key[0]]);
                if offset_reg == None {offset = seq.find(&utils::revcomp(&int_to_minimizer[&key[0]])).unwrap();}
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
        //println!("{}", consensus_str);
        for (kmer, new_str) in bucket_gapped_seqs {
            let mut start_pos = 0;
            for i in 0..new_str.len() {
                if new_str.chars().nth(i).unwrap().to_string() != "-" {
                    start_pos = i;
                    break;
                }


            }
            let mod_str = &consensus_str[start_pos..new_str.len()];
            //println!("{}", mod_str);
            seq_switch.insert(kmer.to_vec(), mod_str.to_string());

        }
        counter += 1;
    }
    seq_switch
}

pub fn buckets_insert(kmer : Vec<u32>, i : usize, buckets : &mut HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>>, mut dbg_nodes: &mut HashMap<Kmer,u32>) {
    for j in 0..kmer.len()-i+1 {
        let mut bucket_idx = Vec::<u32>::new();
        for k in 0..i {
            bucket_idx.push(kmer.to_vec()[j+k]);
        }
        let mut entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<(Vec<u32>, u32)>::new());
        entry.push((kmer.to_vec(), dbg_nodes[&Kmer::make_from(&kmer.to_vec())]));
    }

}