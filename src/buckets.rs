use std::collections::HashMap;
use super::Kmer;
use super::kmer_vec;
use super::Params;

pub fn get_kmers(dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> Vec<Vec<u32>> {
    dbg_nodes.iter().map(|(kmer, count)| kmer_vec::get(kmer).to_vec()).collect()
}

pub fn enumerate_buckets(mut dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params, kmer_seqs: &mut HashMap<Kmer,String>, minim_shift : &mut HashMap<Kmer,(u32,u32)>) -> HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>> {
    let kmers = get_kmers(&mut dbg_nodes, &params);
    let mut buckets : HashMap<Vec<u32>, Vec<(Vec<u32>, u32)>> = HashMap::new();
    for kmer in kmers.iter() {
        buckets_insert(kmer.to_vec(), (params.k-1)/2, &mut buckets, &mut dbg_nodes);
    }
    for (key, entry) in &mut buckets {
        entry.sort_by_key(|tuple| tuple.0.iter().position(|&x| x == key[0]).unwrap());
        entry.reverse();

    }
    buckets
}
pub fn get_consensus(mut dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params, kmer_seqs: &mut HashMap<Kmer,String>, minim_shift : &mut HashMap<Kmer,(u32,u32)>) {
    let mut buckets = enumerate_buckets(&mut dbg_nodes, &params, kmer_seqs, minim_shift);
    println!("{:?}", buckets.len());
    buckets.retain(|key, entry| entry.len() >= 3 && entry.iter().filter(|tuple| tuple.1 == 1).count() == entry.len());
    println!("{:?}", buckets.len());
    for (key, entry) in &buckets {
        let mut max_index = -1;
        let mut min_index = params.k+1;
        let mut entry_collect = Vec::<u32>::new();
        let mut counts : HashMap<u32, u32> = HashMap::new();
        for (kmer, count) in entry.iter() {
            let kmer_node = Kmer::make_from(&kmer.to_vec());
            println!("{:?}\t{}\t{:?}\t{:?}", key, dbg_nodes[&kmer_node], kmer, kmer_seqs[&kmer_node]);
        }
        for (kmer, count) in entry.iter() {
            for lmer in kmer {
                if !entry_collect.contains(lmer) {entry_collect.push(*lmer);}
                let mut count = counts.entry(*lmer).or_insert(0);
                *count += 1;
            }
        }
        let consensus : Vec<u32> = entry_collect.iter().filter(|x| counts[x] > 1).cloned().collect();
        if consensus.len() > params.k {
            let mut multiple = Vec::<Vec<u32>>::new();
            for i in 0..consensus.len()-params.k+1 {
                multiple.push(consensus[i..i+params.k].to_vec());
            }
            for correct in multiple.iter() {
                if correct.len() < params.k {continue;}
                let node = Kmer::make_from(&correct.to_vec());
                let ent = dbg_nodes.entry(node.clone()).or_insert(0);
                *ent += entry.len() as u32 - multiple.len() as u32;
                println!("{:?}", entry_collect);
                println!("{:?}", correct);
                println!("{:?}\t{}\t{:?}\t{:?}\n", key, dbg_nodes[&node], correct, kmer_seqs[&node]);

            }
        }
        else {
            if consensus.len() < params.k {continue;}
            let node = Kmer::make_from(&consensus.to_vec());
            let ent = dbg_nodes.entry(node.clone()).or_insert(0);
            *ent += entry.len() as u32 -1;
            println!("{:?}", entry_collect);
            println!("{:?}", consensus);
            println!("{:?}\t{}\t{:?}\t{:?}\n", key, dbg_nodes[&node], consensus, kmer_seqs[&node]);

        }

    }
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