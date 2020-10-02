use super::Params;
use super::minimizers;
use nthash::ntc64;
use std::collections::HashMap;
use std::collections::VecDeque;
use std::fs::File;
use super::Kmer;
use super::revcomp_aware;
use super::ec_reads;
use super::utils;
use super::poa;
use std::io::BufWriter;
type Buckets<'a> = HashMap<Vec<u64>, Vec<String>>;
#[derive(Clone)]
pub struct Read {
    pub id: String,
    pub minimizers : Vec<String>,
    pub minimizers_pos: Vec<usize>,
    pub transformed: Vec<u64>,
    pub seq: String,
    pub corrected: bool,
}
#[derive(Clone)]
pub struct Lmer {
    pub pos: usize,
    pub hash: u64,
}

impl Read {

    pub fn extract(inp_id: String, inp_seq: String, params: &Params, int_to_minimizer: &HashMap<u64, String>) -> Self {
        if params.w != 0 {
            return Read::extract_windowed(inp_id, inp_seq, params, int_to_minimizer)
        }
        else {
            return Read::extract_density(inp_id, inp_seq, params, int_to_minimizer)
        }
    }

    pub fn extract_windowed(inp_id: String, inp_seq: String, params: &Params, int_to_minimizer: &HashMap<u64, String>) -> Self {
        let l = params.l;
        let w = params.w;
        let mut read_minimizers = Vec::<String>::new();
        let mut read_minimizers_pos = Vec::<usize>::new();
        let mut read_transformed = Vec::<u64>::new();
        let mut Q: VecDeque<Lmer> = VecDeque::new();
        let mut M: VecDeque<Lmer> = VecDeque::new();
        for i in 0..inp_seq.len()-l+1 {
            let mut lmer = &inp_seq[i..i+l];
            let mut lmer = minimizers::normalize_minimizer(&lmer.to_string());
            if lmer.contains("N") {continue;}
            let lmer_hash = ntc64(lmer.as_bytes(), 0, l);
            let lmer_obj = Lmer {pos: i, hash: lmer_hash};
            while !Q.is_empty() && Q.back().unwrap().hash > lmer_hash {
                Q.pop_back();
            }
            Q.push_back(lmer_obj);
            if Q.front().unwrap().pos as i32 <= (i as i32) - (w as i32) {
                while Q.front().unwrap().pos as i32 <= (i as i32) - (w as i32) {
                    Q.pop_front();
                }
                match Q.get(1).cloned() {
                    None => {
                        continue
                    }
                    Some(next) => {
                        while Q.front().unwrap().hash == next.hash {
                            Q.pop_front();
                            if Q.is_empty() {break;}
                        }
                    }
                }
            }
            match M.back() {
                None => {
                    M.push_back(Q.front().unwrap().clone());
                }
                Some(back) => {
                    if Q.is_empty() {continue;}
                    if back.hash != Q.front().unwrap().hash {
                        M.push_back(Q.front().unwrap().clone());
                    }
                }
            }
           
            
        }
        let mut read_minimizers : Vec::<String> = M.iter().map(|lmer| int_to_minimizer[&lmer.hash].to_string()).collect();
        let mut read_minimizers_pos : Vec::<usize> = M.iter().map(|lmer| lmer.pos).collect();
        let mut read_transformed : Vec::<u64> = M.iter().map(|lmer| lmer.hash).collect();

        Read {id: inp_id, minimizers: read_minimizers, minimizers_pos: read_minimizers_pos, transformed: read_transformed, seq: inp_seq, corrected: false}

    }
    pub fn extract_density(inp_id: String, inp_seq: String, params: &Params, int_to_minimizer : &HashMap<u64, String>) -> Self {
        let size_miniverse = params.size_miniverse as u64;
        let density = params.density;
        let l = params.l;
        let mut read_minimizers = Vec::<String>::new();
        let mut read_minimizers_pos = Vec::<usize>::new();
        let mut read_transformed = Vec::<u64>::new();
        for i in 0..inp_seq.len()-l+1 {
            let mut lmer = &inp_seq[i..i+l];
            let mut lmer = minimizers::normalize_minimizer(&lmer.to_string());
            let hash = ntc64(lmer.as_bytes(), 0, l);
            if (!lmer.contains("N")) && (hash as f64) < u64::max_value() as f64 * density/(l as f64) {
                read_minimizers.push(lmer.to_string());
                read_minimizers_pos.push(i);
                read_transformed.push(hash);
            }
        }
        Read {id: inp_id, minimizers: read_minimizers, minimizers_pos: read_minimizers_pos, transformed: read_transformed, seq: inp_seq, corrected: false}
    }

    pub fn label(&self, read_seq: String, read_minimizers: Vec<String>, read_minimizers_pos: Vec<usize>, read_transformed: Vec<u64>, corrected_map: &mut HashMap<String, (String, Vec<String>, Vec<usize>, Vec<u64>)>) {
        corrected_map.insert(self.id.to_string(), (read_seq, read_minimizers, read_minimizers_pos, read_transformed));
    }

    pub fn read_to_kmers(&mut self, kmer_origin: &mut HashMap<Kmer,String>, dbg_nodes: &mut HashMap<Kmer,u32> , kmer_seqs: &mut HashMap<Kmer,String>, minim_shift : &mut HashMap<Kmer, (usize, usize)>, params: &Params) {
        let k = params.k;
        let l = params.l;
        let n = params.n;
        let min_kmer_abundance = params.min_kmer_abundance;
        let levenshtein_minimizers = params.levenshtein_minimizers;
        let output_base_space = params.output_base_space;
        for i in 0..(self.transformed.len()-k+1) {
            let mut node : Kmer = Kmer::make_from(&self.transformed[i..i+k]);
            let mut seq_reversed = false;
            if revcomp_aware { 
                let (node_norm, reversed) = node.normalize(); 
                node = node_norm;
                seq_reversed = reversed;
            } 
            let entry = dbg_nodes.entry(node.clone()).or_insert(0);
            *entry += 1;

            if ! output_base_space { continue;}


            //if *entry == min_kmer_abundance as u32 {
                let mut seq = self.seq[self.minimizers_pos[i] as usize..(self.minimizers_pos[i+k-1] as usize + l)].to_string();
                if seq_reversed {
                    seq = utils::revcomp(&seq);
                }
                kmer_seqs.insert(node.clone(), seq.clone());
                let origin = format!("{}_{}_{}", self.id, self.minimizers_pos[i].to_string(), self.minimizers_pos[i+k-1].to_string());
                kmer_origin.insert(node.clone(), origin);
                let position_of_second_minimizer = match seq_reversed {
                    true => self.minimizers_pos[i+k-1]-self.minimizers_pos[i+k-2],
                    false => self.minimizers_pos[i+1]-self.minimizers_pos[i]
                };
                let position_of_second_to_last_minimizer = match seq_reversed {
                    true => self.minimizers_pos[i+1]-self.minimizers_pos[i],
                    false => self.minimizers_pos[i+k-1]-self.minimizers_pos[i+k-2]
                };
                minim_shift.insert(node.clone(), (position_of_second_minimizer, position_of_second_to_last_minimizer));
                if levenshtein_minimizers == 0 {
                    for minim in &self.minimizers[i..i+k-1] {
                        debug_assert!((!&seq.find(minim).is_none()) || (!utils::revcomp(&seq).find(minim).is_none()));
                    }
                }
            //}
        }
    }

    
    //pub fn write_to_poa
    pub fn query(&mut self, int_to_minimizer: &HashMap<u64, String>, poa_map: &mut HashMap<String, Vec<String>>, buckets: &Buckets, params : &Params, mut corrected_map: &mut HashMap<String, (String, Vec<String>, Vec<usize>, Vec<u64>)>, reads_by_id: &HashMap<String, Read>) {
        let n = params.n;
        let k = params.k;
        let l = params.l;
        let mut read_minimizers = &self.minimizers;
        let mut read_minimizers_pos = &self.minimizers_pos;
        let mut read_transformed = &self.transformed;
        let mut seq_id = &self.id;
        let mut seq_str = &self.seq;
        let mut scoring = poa::Scoring::new(-1, 0, |a: u64, b: u64| if a == b { 1i32 } else { -1i32 });
        let mut aligner = poa::Aligner::new(scoring, &read_transformed);
        let mut aligned : HashMap<&Vec<u64>, bool> = HashMap::new();
        let mut bucket_reads = Vec::<&Read>::new();
        let mut poa_ids = Vec::<String>::new();
        let mut pair_map : HashMap<(u64, u64), String> = HashMap::new();
        for i in 0..read_transformed.len()-1 {
            pair_map.insert((read_transformed[i], read_transformed[i+1]), seq_str[read_minimizers_pos[i] as usize..read_minimizers_pos[i+1] as usize].to_string());
        }
        for i in 0..read_transformed.len()-n+1 {
            let bucket_idx = &read_transformed[i..i+n];
            let entry = &buckets[bucket_idx];
            for id in entry.iter() {
                let query = &reads_by_id[id];
                if &query.transformed == read_transformed {continue;}
                let aligned_entry = aligned.entry(&query.transformed).or_insert(false);
                if !*aligned_entry {
                    bucket_reads.push(query);
                    *aligned_entry = true;  
                }
            }   
        }
        let mut bucket_reads : Vec<(&Read, f64)> = bucket_reads.iter().map(|seq| (*seq, minimizers::dist(self, seq, &params))).filter(|(seq, dist)| *dist < 0.25).collect();
        bucket_reads.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        for (read, dist) in bucket_reads.iter() {
            let transformed = &read.transformed;
            let seq = &read.seq;
            let pos = &read.minimizers_pos;
            for i in 0..transformed.len()-1 {
                pair_map.insert((transformed[i], transformed[i+1]), seq[pos[i] as usize ..pos[i+1] as usize].to_string());
            }
        }
        let mut max_len = 100;
        if bucket_reads.len() > max_len {bucket_reads = bucket_reads[0..max_len].to_vec();}
        //if bucket_reads.len() == 0 {println!("Read has no neighbors");}
        for i in 0..bucket_reads.len() {
            poa_ids.push(bucket_reads[i].0.id.to_string());
            aligner.global(&bucket_reads[i].0.transformed);
            aligner.add_to_graph();
        }
        let mut consensus = aligner.poa.consensus(&params);
        let consensus_read = consensus.iter().map(|minim| int_to_minimizer[minim].to_string()).collect::<Vec<String>>();
        let mut consensus_str = String::new();
        let mut pos_idx = 0;
        let mut consensus_pos = Vec::<usize>::new();
        if consensus.len() == 0 {return}
        for i in 0..consensus.len()-1 {
            consensus_pos.push(pos_idx);
            let mut insert = String::new();
            if !pair_map.contains_key(&(consensus[i], consensus[i+1])) {
                let alts : Vec<String> = pair_map.clone().into_iter().filter(|(a, b)| a.0 == consensus[i]).map(|tuple| tuple.1).collect();
                insert = int_to_minimizer[&consensus[i]].to_string();
                for _ in 0..alts[0].len() {insert.push_str("N")};
            }
            else { insert = pair_map[&(consensus[i], consensus[i+1])].to_string();}
            consensus_str.push_str(&insert);
            pos_idx += insert.len();
        }
        consensus_pos.push(pos_idx);
        consensus_str.push_str(&int_to_minimizer[&consensus[consensus.len()-1]]);
        let mut corrected_count = 0;
        let mut threshold = params.correction_threshold;
        if params.correction_threshold == 0 {threshold = 0 as i32;}
        for (read, dist) in bucket_reads.iter_mut() {
            if corrected_count >= threshold {break;}
            if !read.corrected {
                read.label(consensus_str.to_string(), consensus_read.to_vec(), consensus_pos.to_vec(), consensus.to_vec(), &mut corrected_map);
                corrected_count += 1;
            }
        }
        poa_map.insert(seq_id.to_string(), poa_ids.to_vec());
        self.seq = consensus_str;
        self.minimizers = consensus_read;
        self.minimizers_pos = consensus_pos;
        self.transformed = consensus;
        self.corrected = true;
    }
}
