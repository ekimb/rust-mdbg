use std::collections::HashMap;
use super::Kmer;
use super::kmer_vec;
use super::Params;

pub fn levenshtein(a: &Vec<u32>, b: &Vec<u32>) -> usize {

    let mut result = 0;

    /* Shortcut optimizations / degenerate cases. */
    if a == b {
        return result;
    }
    
    let length_a = a.len();
    let length_b = b.len();

    if length_a == 0 {
        return length_b;
    }

    if length_b == 0 {
        return length_a;
    }

    /* Initialize the vector.
    *
    * This is why it’s fast, normally a matrix is used,
    * here we use a single vector. */
    let mut cache: Vec<usize> = vec![0; length_a];
    let mut index_a = 0;
    let mut distance_a;
    let mut distance_b;

    while index_a < length_a {
        index_a += 1;
        cache[index_a - 1] = index_a;
    }

    /* Loop. */
    for (index_b, code_b) in b.iter().enumerate() {
        result = index_b;
        distance_a = index_b;

        for (index_a, code_a) in a.iter().enumerate() {
            distance_b = if code_a == code_b {
                distance_a
            } else {
                distance_a + 1
            };

            distance_a = cache[index_a];

            result = if distance_a > result {
                if distance_b > result {
                    result + 1
                } else {
                    distance_b
                }
            } else {
                if distance_b > distance_a {
                    distance_a + 1
                } else {
                    distance_b
                }
            };

            cache[index_a] = result;
        }
    }

    result
}


pub fn get_nonsolid(dbg_nodes: &mut HashMap<Kmer,u32>, params : &Params) -> Vec<Vec<u32>> {
    dbg_nodes.iter().map(|(kmer, count)| kmer_vec::get(kmer).to_vec()).collect()
}
pub fn dbscan(kmers : &Vec<Vec<u32>>, epsilon : u32) ->  HashMap<Vec::<u32>,Vec<Vec<u32>>> {
    let mut neighbors_map : HashMap<Vec::<u32>,Vec<Vec<u32>>> = HashMap::new();
    let mut label : HashMap<Vec::<u32>,String> = HashMap::new();
    let mut cluster_count = 0;
    for kmer in kmers.iter() {
        if label.contains_key(kmer) {continue;}
        let mut neighbors : Vec::<Vec<u32>> = range_query(kmers, kmer, epsilon);
        if neighbors.len() < 3 {
            label.insert(kmer.to_vec(), "N".to_string());
            //println!("Inserted noise");
            continue;
        }
        cluster_count += 1;
        label.insert(kmer.to_vec(), "C".to_string());
        //println!("Inserted cluster");
        let mut seeds : Vec::<Vec<u32>> = neighbors.iter().filter(|x| x != &kmer).map(|x| x.to_vec()).collect();
        for seed in seeds.to_vec().iter() {
            if label.contains_key(seed) && label[seed] == "N".to_string() {label.insert(seed.to_vec(), "C".to_string());}
            if label.contains_key(seed) {continue;}
            label.insert(seed.to_vec(), "C".to_string());
            let mut seed_neighbors = range_query(kmers, seed, epsilon);
            if seed_neighbors.len() >= 3 {
                seeds.extend(neighbors.to_vec());
            }
        }
        neighbors_map.insert(kmer.to_vec(), neighbors);
    }
    neighbors_map
}
pub fn range_query(kmers : &Vec<Vec<u32>>, kmer : &Vec<u32>, epsilon : u32) -> Vec<Vec<u32>> {
    let mut neighbors = Vec::<Vec<u32>>::new();
    for other_kmer in kmers.iter() {
        if levenshtein(kmer, other_kmer) <= epsilon as usize && other_kmer != kmer {
            neighbors.push(other_kmer.to_vec())
        }
    }
    neighbors
}

pub fn cluster_kmers(mut dbg_nodes: &mut HashMap<Kmer,u32> , kmer_seqs: &mut HashMap<Kmer,String> , minim_shift : &mut HashMap<Kmer,(u32,u32)>, params : &Params) {
    let kmers = get_nonsolid(&mut dbg_nodes, &params);
    println!("Got nonsolid kmers");
    let mut neighbors_map = dbscan(&kmers, 1);
    for (kmer, neighbors) in &neighbors_map {
        let mut consensus = Vec::<u32>::new();
        for i in 0..params.k {
            let mut counts : HashMap<u32, u32> = HashMap::new();
            let mut lmers = Vec::<u32>::new();
            for neighbor in neighbors.iter() {
                lmers.push(neighbor[i]);
                let solid_kmer = Kmer::make_from(&neighbor);
                if dbg_nodes.contains_key(&solid_kmer) && dbg_nodes[&solid_kmer] >= 2 as u32 {
                    for j in 0..dbg_nodes[&solid_kmer] {lmers.push(neighbor[i]);}
                }
            }
            for lmer in lmers.iter() {
                let mut count = counts.entry(*lmer).or_insert(0);
                *count += 1;
            }
            let mut max = counts.clone().into_iter().max_by_key(|&(tuple)| tuple.1).unwrap();
            consensus.push(max.0);
        }
        let mut consensus_node : Kmer = Kmer::make_from(&consensus);
        let (node_norm, reversed) = consensus_node.normalize(); 
        consensus_node = node_norm;
        let seq_reversed = reversed;
        println!("Kmer {:?}, neighbors {:?}", kmer, neighbors);
        if dbg_nodes.contains_key(&consensus_node) && dbg_nodes[&consensus_node] >= 2 as u32 {
            println!("Found existing kmer seq {:?}", kmer_seqs[&consensus_node]);
            for neighbor in neighbors {
                let mut neighbor_node : Kmer = Kmer::make_from(&neighbor);
                let (node_norm, reversed) = neighbor_node.normalize(); 
                neighbor_node = node_norm;
                let seq_reversed = reversed;
                let mut count = dbg_nodes.entry(consensus_node.clone()).or_insert(0);
                *count += 1;
                if dbg_nodes.contains_key(&neighbor_node) && dbg_nodes[&neighbor_node] < 2 as u32 {
                    dbg_nodes.remove(&neighbor_node);
                    println!("Dbg removed");
                    if kmer_seqs.contains_key(&neighbor_node) {
                        kmer_seqs.remove(&neighbor_node);
                        println!("Kmerseq removed");

                    }
                    if minim_shift.contains_key(&neighbor_node){
                        minim_shift.remove(&neighbor_node);
                        println!("Shift removed");

                    }
                    
                }
                
            }

        }
        else {
            println!("Found new kmer {:?}", consensus)
        }
    }
}


