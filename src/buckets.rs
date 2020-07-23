use std::collections::HashMap;
use super::Kmer;
use super::kmer_vec;
use super::Params;
use super::utils;
use super::extract_minimizers;

pub fn get_kmers(dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> Vec<Vec<u32>> {
    dbg_nodes.iter().map(|(kmer, count)| kmer_vec::get(kmer).to_vec()).collect()
}

pub fn enumerate_buckets(mut seq_mins : &mut Vec<Vec<u32>>, mut dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> HashMap<Vec<u32>, Vec<Vec<u32>>> {
    let mut buckets : HashMap<Vec<u32>, Vec<Vec<u32>>> = HashMap::new();
    for seq in seq_mins.iter() {
            buckets_insert(seq.to_vec(), 2, &mut buckets);
    }
    for (key, entry) in &mut buckets {
        entry.sort_by_key(|seq| seq.iter().position(|&x| x == key[0]).unwrap());
        entry.reverse();

    }
    println!("{:?}", buckets.len());
    //buckets.retain(|key, entry| entry.len());
    println!("{:?}", buckets.len());
    buckets
}

pub fn query_buckets(read_transformed : Vec<u32>, mut dbg_nodes: &mut HashMap<Kmer,u32>, buckets : &mut HashMap<Vec<u32>, Vec<Vec<u32>>>, mut seq_str: &mut String, params : &Params, mut seq_mins: &mut Vec<Vec<u32>>, lmer_counts: &HashMap<String,u32>, minimizer_to_int : &HashMap<String,u32>, int_to_minimizer : &HashMap<u32,String>) -> Vec<u32>{
    let sub = 2;
    let k = params.k;
    let mut consensus_seqs = Vec::<Vec<u32>>::new();
    let mut potential_kmers : HashMap<u32, Vec<Vec<u32>>> = HashMap::new();
    //print!("\n");
    for i in 0..read_transformed.len()-k+1 {
        let mut candidate_kmers = Vec::<Vec<u32>>::new();
        let mut og_kmer = read_transformed[i..i+k].to_vec();
        //println!("Kmer {:?}", og_kmer);
        for j in 0..k-sub+1 {
            let mut bucket_idx = Vec::<u32>::new();
            bucket_idx = og_kmer[j..j+sub].to_vec();
            let mut entry = Vec::<Vec<u32>>::new();
            if buckets.contains_key(&bucket_idx) {
                entry = buckets[&bucket_idx].to_vec();
            }
           // println!("Here with {:?}", bucket_idx);
            let mut bucket_idx_rev = bucket_idx.to_vec().clone();
            bucket_idx_rev.reverse();
            let mut entry_rev = Vec::<Vec<u32>>::new();
            if buckets.contains_key(&bucket_idx_rev) {
                entry_rev = buckets[&bucket_idx_rev].to_vec();
                //println!("Found revcomp");
            } 
            let mut bucket_seqs = Vec::<Vec<u32>>::new();
            if entry.len() != 0 {
                for seq_min in entry.iter() {
                    bucket_seqs.push(seq_min.to_vec());
                }
            }
            if entry_rev.len() != 0 {
                for seq_min in entry_rev.iter() {
                    bucket_seqs.push(seq_min.to_vec());
                }
            }
            let mut bucket_gapped_seqs = Vec::<Vec<u32>>::new();
            let mut counts : HashMap<u32, u32> = HashMap::new();
            let mut min_counts : HashMap<u32, Vec<u32>> = HashMap::new();
            // println!("Found false kmer");
            let mut prev_len = 0;
            let mut min_prev_len = 99999;
            let mut og_str_idx = (0, false);
            for seq in bucket_seqs.iter() {
                let mut new_seq = Vec::<u32>::new();
                if seq.to_vec() == og_kmer {og_str_idx = (bucket_seqs.iter().position(|x| x.to_vec() == og_kmer).unwrap(), false)};
                let mut offset;
                let mut offset_reg = seq.iter().position(|&x| x == bucket_idx[0]);
                offset = offset_reg.unwrap();
                if offset > prev_len {prev_len = offset;}
                if offset < min_prev_len {min_prev_len = offset;}       
            }
            for seq in bucket_seqs.iter() {
                //print!("Seq\t");
                //for min in seq.iter() {
                //    print!("{}\t", min);
            // }
                //print!("\n");
                    let mut new_seq = Vec::<u32>::new();
                    let mut offset;
                    let mut offset_reg = seq.iter().position(|&x| x == bucket_idx[0]);
                    offset = offset_reg.unwrap();
                let mut new_seq = Vec::<u32>::new();
                if prev_len != offset {
                    for i in offset..prev_len {
                        new_seq.push(0)
                    }
                }       
                for min in seq.iter() {
                    new_seq.push(*min);
                }     
               // print!("New\t");
                for min in new_seq.iter() {
                 //   print!("{}\t", min);
                }
                //print!("\n");
                bucket_gapped_seqs.push(new_seq.to_vec());
            }
            //println!("\n");
            for seq in bucket_gapped_seqs.iter() {
                for j in 0..seq.len() {
                    let mut pos_entry = min_counts.entry(j as u32).or_insert(Vec::<u32>::new());
                    pos_entry.push(seq[j]);
                }     
            }
            let mut consensus_seq = Vec::<u32>::new();
            for pos in 0..min_counts.keys().len() {
                let mut min_freq : HashMap<u32, u32> = HashMap::new();
                for min in min_counts[&(pos as u32)].iter() {
                    if *min == 0 {continue;}
                    let min_cnt = min_freq.entry(*min).or_insert(0);
                    *min_cnt += 1;
                }
                let max_min = min_freq.into_iter().max_by_key(|&(_, count)| count).map(|(val, _)| val).unwrap();
                consensus_seq.push(max_min);
            }
            //println!("Consensus seq {}", consensus_seq.len());
            for min in consensus_seq.iter() {
                    //print!("{}\t", min);
                }
               // print!("\n");
            let mut seq_min_new = Vec::<u32>::new();
            for seq in bucket_gapped_seqs.iter() {
                if bucket_gapped_seqs.iter().position(|x| x.to_vec() == seq.to_vec()).unwrap() == og_str_idx.0 {
                    //println!("gapped seq is {}", og_str_idx.0);
                    let mut start_pos = 0;
                    for i in 0..seq.len() {
                        if seq[i] != 0 {
                            start_pos = i;
                            break;
                        }
                    }
                    seq_min_new = consensus_seq[start_pos..seq.len()].to_vec();
                    if !candidate_kmers.contains(&seq_min_new.to_vec()) {
                        candidate_kmers.push(seq_min_new.to_vec());
                        //println!("{:?}", read_transformed);
                        //println!("Now");
                        //println!("{:?}", seq_min_new);
                        break;
                    }
                }
            }
        }
       // println!("Candidate kmers for pos {}", i);
        for kmer in candidate_kmers.iter() {
            let pos_entry = potential_kmers.entry(i as u32).or_insert(Vec::<Vec<u32>>::new());
            pos_entry.push(kmer.to_vec());
        for min in kmer.iter() {
            //print!("{}\t", min);
            }
            //print!("\n");
        }
    }
    for pos in 0..potential_kmers.keys().len() {
        println!("Kmers for pos {}", pos);
        for kmer in potential_kmers[&(pos as u32)].iter() {
            println!("{:?}", kmer);
        }
    }
    find_path(&read_transformed, &potential_kmers, &params);
    if consensus_seqs.len() > 2 {
        vec![1]
        
    }
    else {
        read_transformed.to_vec()
    }

}

pub fn recur_path(mins_len : usize, potential_kmers : &HashMap<u32, Vec<Vec<u32>>>, mut build : &mut Vec<u32>, pos : u32, params : &Params) -> Vec<u32> {
    let k = params.k;
    let pos_new = pos + 1;
    let mut extra_kmer = 0;
    if potential_kmers.contains_key(&pos_new) {
        for kmer in potential_kmers[&pos_new].iter() {
            let build_suffix = build[build.len()-k+1..build.len()].to_vec();
            let kmer_prefix = kmer[0..k-1].to_vec();
            if build_suffix == kmer_prefix {
                println!("Pos new {}", pos_new);
                build.push(kmer[k-1]);
                println!("build {:?}", build);
                *build = recur_path(mins_len, &potential_kmers, &mut build, pos_new, &params);                
            }
        }
    }
    else {
        return vec![1];
    }
    return build.to_vec()
}
pub fn find_path(read_transformed : &Vec<u32>, potential_kmers : &HashMap<u32, Vec<Vec<u32>>>, params : &Params) {
    let mut build = Vec::<u32>::new();
    let mins_len = read_transformed.len();
    for kmer in potential_kmers[&0].iter() {
            build = kmer.clone();
            build = recur_path(mins_len, &potential_kmers, &mut build, 0, &params).to_vec();
            println!("Build {:?}", build);

    }


}
pub fn buckets_insert(seq : Vec<u32>, i : usize, buckets : &mut HashMap<Vec<u32>, Vec<Vec<u32>>>) {
    for j in 0..seq.len()-i+1 {
        let mut bucket_idx = Vec::<u32>::new();
        for k in 0..i {
            bucket_idx.push(seq.to_vec()[j+k]);
        }
        let mut entry = buckets.entry(bucket_idx.to_vec()).or_insert(Vec::<Vec<u32>>::new());
        entry.push((seq.to_vec()));
    }

}