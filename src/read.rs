use super::{Kmer, Repr, CorrMap, Params, minimizers, RacyBloom, REVCOMP_AWARE, utils, poa};
use nthash::{NtHashIterator};
use std::collections::{HashMap,HashSet};
//use std::collections::hash_map::Entry::{Occupied, Vacant};
//use std::collections::VecDeque;
use super::utils::pretty_minvec;
use std::collections::VecDeque;

type Buckets<'a> = HashMap<Vec<u64>, Vec<String>>;

#[derive(Clone, Default)]
pub struct Read {
    pub id: String,
    pub minimizers : Vec<String>,
    pub minimizers_pos: Vec<usize>,
    pub transformed: Vec<u64>,
    pub seq: String, 
    //pub seq_str: &'a str, // an attempt to avoid copying the string sequence returned by the fasta parser (seems too much effort to implement for now)
    pub corrected: bool
}
#[derive(Clone)]
pub struct Lmer {
    pub pos: usize,
    pub hash: u64
}


const SEQ_NT4_TABLE: [u8; 256] =
   [0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4];

// copy from http://www.cse.yorku.ca/~oz/hash.html:

pub fn hash(mut key: u64, mask: u64) -> u64 {
        key = (!key + (key << 21)) & mask;
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8)) & mask;
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4)) & mask;
        key = key ^ key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
}


pub fn update_window(q: &mut VecDeque<u64>, q_pos: &mut VecDeque<usize>, q_min_val: u64, q_min_pos: i32, new_strobe_hashval: u64, i: usize, new_minimizer: bool) -> (u64, i32, bool) {
    q.pop_front();
    let popped_index = q_pos.pop_front();
    q.push_back(new_strobe_hashval);
    q_pos.push_back(i);
    let mut min_val = q_min_val;
    let mut min_pos = q_min_pos;
    let mut new_minim = new_minimizer;
    if min_pos == popped_index.unwrap() as i32 {
        min_val = u64::max_value();
        min_pos = i as i32;
        for j in (0..q.len()).rev() {
            if q[j] < min_val {
                min_val = q[j];
                min_pos = q_pos[j] as i32;
                new_minim = true;
            }
        }
    }
    else if new_strobe_hashval < min_val { // the new value added to queue is the new minimum
        min_val = new_strobe_hashval;
        min_pos = i as i32;
        new_minim = true;
    }
    (min_val, min_pos, new_minim)
}



impl Read {
    pub fn extract(inp_id: &str, inp_seq: String, params: &Params, minimizer_to_int: &HashMap<String, u64>, uhs_bloom: &RacyBloom, lcp_bloom: &RacyBloom) -> Self {
        if params.uhs {Read::extract_uhs(inp_id, inp_seq, params, minimizer_to_int, uhs_bloom)}
        else if params.lcp {Read::extract_lcp(inp_id, inp_seq, params, minimizer_to_int, lcp_bloom)}
        else if params.use_syncmers { Read::extract_syncmers(inp_id, inp_seq, params) }
	else {Read::extract_density(inp_id, inp_seq, params, minimizer_to_int)}
    }

    //delete
    pub fn extract_lcp(inp_id: &str, inp_seq_raw: String, params: &Params, minimizer_to_int: &HashMap<String, u64>, lcp_bloom: &RacyBloom) -> Self {
        let density = params.density;
        let l = params.l;
        let read_minimizers = Vec::<String>::new();
        let mut read_minimizers_pos = Vec::<usize>::new();
        let mut read_transformed = Vec::<u64>::new();
        let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;
        let tup;
        let inp_seq;
        if !params.use_hpc {
            tup = Read::encode_rle(&inp_seq_raw); //get HPC sequence and positions in the raw nonHPCd sequence
            inp_seq = tup.0; //assign new HPCd sequence as input
        }
        else {
            inp_seq = inp_seq_raw.clone(); //already HPCd before so get the raw sequence
        }
        let iter = NtHashIterator::new(inp_seq.as_bytes(), l).unwrap().enumerate().filter(|(_i, x)| *x <= hash_bound);
        for (i, mut hash) in iter {
            let lmer = &inp_seq[i..i+l];
            if lmer.contains('N') {continue;}
            if params.error_correct || params.has_lmer_counts {
                let res = minimizer_to_int.get(lmer); // allows to take the 'skip' array into account
                if res.is_none() {continue;} // possible discrepancy between what's calculated in minimizers_preparation() and here
                hash = *res.unwrap();
            }
            if lcp_bloom.get().check_and_add(&hash) {
                read_minimizers_pos.push(i);
                read_transformed.push(hash);
            }
        }
        Read {id: inp_id.to_string(), minimizers: read_minimizers, minimizers_pos: read_minimizers_pos, transformed: read_transformed, seq: inp_seq_raw, corrected: false}
    }
    pub fn extract_uhs(inp_id: &str, inp_seq_raw: String, params: &Params, minimizer_to_int: &HashMap<String, u64>, uhs_bloom: &RacyBloom) -> Self {
        let density = params.density;
        let l = params.l;
        let read_minimizers = Vec::<String>::new();
        let mut read_minimizers_pos = Vec::<usize>::new();
        let mut read_transformed = Vec::<u64>::new();
        let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;
        let tup;
        let inp_seq;
        if !params.use_hpc {
            tup = Read::encode_rle(&inp_seq_raw); //get HPC sequence and positions in the raw nonHPCd sequence
            inp_seq = tup.0; //assign new HPCd sequence as input
        }
        else {
            inp_seq = inp_seq_raw.clone(); //already HPCd before so get the raw sequence
        }
        let iter = NtHashIterator::new(inp_seq.as_bytes(), l).unwrap().enumerate().filter(|(_i, x)| *x <= hash_bound);
        for (i, hash) in iter {
            let lmer = &inp_seq[i..i+l];
            let mut hash : u64 = hash;
            if params.error_correct || params.has_lmer_counts {
                let res = minimizer_to_int.get(lmer); // allows to take the 'skip' array into account
                if res.is_none() {continue;} // possible discrepancy between what's calculated in minimizers_preparation() and here
                hash = *res.unwrap();
            }
            if uhs_bloom.get().check_and_add(&hash) {
                read_minimizers_pos.push(i);
                read_transformed.push(hash);
            }
        }
        Read {id: inp_id.to_string(), minimizers: read_minimizers, minimizers_pos: read_minimizers_pos, transformed: read_transformed, seq: inp_seq_raw, corrected: false}
    }
    pub fn encode_rle(inp_seq: &str) -> (String, Vec<usize>) {
        let mut prev_char = '#';
        let mut hpc_seq = String::new();
        let mut pos_vec = Vec::<usize>::new();
        let mut prev_i = 0;
        for (i, c) in inp_seq.chars().enumerate() {
            if c == prev_char && "ACTGactgNn".contains(c) {continue;}
            if prev_char != '#' {
                hpc_seq.push(prev_char);
                pos_vec.push(prev_i);
                prev_i = i;
            }
            prev_char = c;
        }
        hpc_seq.push(prev_char);
        pos_vec.push(prev_i);
        (hpc_seq, pos_vec)
    }
    
    pub fn extract_density(inp_id: &str, inp_seq_raw: String, params: &Params, minimizer_to_int: &HashMap<String, u64>) -> Self {
        let density = params.density;
        let l = params.l;
        let read_minimizers = Vec::<String>::new();
        let mut read_minimizers_pos = Vec::<usize>::new();
        let mut read_transformed = Vec::<u64>::new();
        //println!("parsing new read: {}\n",inp_seq);
        let hash_bound = ((density as f64) * (u64::max_value() as f64)) as u64;
        let mut tup = (String::new(), Vec::<usize>::new());
        let inp_seq;
        if !params.use_hpc {
            tup = Read::encode_rle(&inp_seq_raw); //get HPC sequence and positions in the raw nonHPCd sequence
            inp_seq = tup.0; //assign new HPCd sequence as input
        }
        else {
            inp_seq = inp_seq_raw.clone(); //already HPCd before so get the raw sequence
        }
        if inp_seq.len() < l {
            return Read {id: inp_id.to_string(), minimizers: read_minimizers, minimizers_pos: read_minimizers_pos, transformed: read_transformed, seq: inp_seq_raw, corrected: false};
        }
        let iter = NtHashIterator::new(inp_seq.as_bytes(), l).unwrap().enumerate().filter(|(_i, x)| *x <= hash_bound);
        for (i,hash) in iter {
            let lmer = &inp_seq[i..i+l];
            let mut hash :u64 = hash;
            if params.error_correct || params.has_lmer_counts {
                let res = minimizer_to_int.get(lmer); // allows to take the 'skip' array into account
                if res.is_none() {continue;} // possible discrepancy between what's calculated in minimizers_preparation() and here
                hash = *res.unwrap();
            }
            if !params.use_hpc {read_minimizers_pos.push(tup.1[i]);} //if not HPCd need raw sequence positions
            else {read_minimizers_pos.push(i);} //already HPCd so positions are the same
            read_transformed.push(hash);
        }
        Read {id: inp_id.to_string(), minimizers: read_minimizers, minimizers_pos: read_minimizers_pos, transformed: read_transformed, seq: inp_seq_raw, corrected: false}
    }

    // copied from hifimap's code
    // except that we use down-sampled syncmers here, otherwise we'd get too many minimizers
    pub fn extract_syncmers(inp_id: &str, inp_seq_raw: String, params: &Params) -> Self {
        let l = params.l;
	let hash_bound = ((params.density as f64) * 4_usize.pow(l as u32) as f64) as u64;


        // boilerplate code for all minimizer schemes
        let mut tup = (String::new(), Vec::<usize>::new());
        let seq;
        if !params.use_hpc {
            tup = Read::encode_rle(&inp_seq_raw); //get HPC sequence and positions in the raw nonHPCd sequence
            seq = tup.0.as_bytes(); //assign new HPCd sequence as input
        }
        else {
            seq = inp_seq_raw.as_bytes(); //already HPCd before so get the raw sequence
        }
        if seq.len() < l {
            let read_minimizers = Vec::<String>::new();
            let read_minimizers_pos = Vec::<usize>::new();
            let read_transformed = Vec::<u64>::new();
            return Read {id: inp_id.to_string(), minimizers: read_minimizers, minimizers_pos: read_minimizers_pos, transformed: read_transformed, seq: inp_seq_raw, corrected: false};
        }

	let s = params.s;
	//let wmin = params.wmin;
	//let wmax = params.wmax;
	// using hifimap defaults
	//let wmin = l/(l-s+1)+2; // actually unused here
	//let wmax = l/(l-s+1)+10;
	let smask : u64 = ((1 as u64) << 2*s) - 1;
	let lmask : u64 = ((1 as u64) << 2*l) - 1;
	let t = f64::ceil((l - s + 1) as f64 / 2.0) as usize;
        let mut seq_hashes = Vec::new();
        let mut pos_to_seq_coord = Vec::new();
        let mut qs = VecDeque::<u64>::new();
        let mut qs_pos = VecDeque::<usize>::new();
        let seq_len = seq.len();
        let mut qs_size = 0;
        let mut qs_min_val = u64::max_value();
        let mut qs_min_pos : i32 = -1;
        let mut xl : [u64; 2] = [0; 2];
        let mut xs : [u64; 2] = [0; 2];
        let mut lp = 0;
        let lshift : u64 = (l as u64 - 1) * 2;
        let sshift : u64 = (s as u64 - 1) * 2;
        for i in 0..seq_len {
            let c = SEQ_NT4_TABLE[seq[i] as usize];
            if c < 4 {
                xl[0] = (xl[0] << 2 | c as u64) & lmask;
                xl[1] = xl[1] >> 2 | ((3 - c) as u64) << lshift;
                xs[0] = (xs[0] << 2 | c as u64) & smask;
                xs[1] = xs[1] >> 2 | ((3 - c) as u64) << sshift;
                lp += 1;
                if s != 0 { //ksyncmer or kstrobemer
                    if lp >= s {
                        let ys : u64 = match xs[0] < xs[1]{
                            true => xs[0],
                            false => xs[1]
                        };
                        let hash_s = hash(ys, smask);
                        if qs_size < l - s {
                            qs.push_back(hash_s);
                            qs_pos.push_back(i - s + 1);
                            qs_size += 1;
                        }
                        else if qs_size == l - s {
                            qs.push_back(hash_s);
                            qs_pos.push_back(i - s + 1);
                            qs_size += 1;
                            for j in 0..qs_size {
                                if qs[j] < qs_min_val {
                                    qs_min_val = qs[j];
                                    qs_min_pos = qs_pos[j] as i32;
                                }
                            }
                            if qs_min_pos == qs_pos[t-1] as i32 {
                                let yl : u64 = match xl[0] < xl[1]{
                                    true => xl[0],
                                    false => xl[1]
                                };
                                let hash_l = hash(yl, lmask);
                                if hash_l <= hash_bound {
                                    seq_hashes.push(hash_l);
                                    //pos_to_seq_coord.push(i - l + 1);
                                    if !params.use_hpc {pos_to_seq_coord.push(tup.1[i-l+1]);} //if not HPCd need raw sequence positions
                                    else {pos_to_seq_coord.push(i-l+1);} //already HPCd so positions are the same
                                }
                            }
                        }
                        else {
                            //let mut new_minimizer = false; // is never used
                            let tuple = update_window(&mut qs, &mut qs_pos, qs_min_val, qs_min_pos, hash_s, i - s + 1, /*new_minimizer*/ false);
                            qs_min_val = tuple.0; qs_min_pos = tuple.1; 
                            // new_minimizer = tuple.2;// is never used
                            if qs_min_pos == qs_pos[t-1] as i32 {
                                let yl : u64 = match xl[0] < xl[1] {
                                    true => xl[0],
                                    false => xl[1]
                                };
                                let hash_l = hash(yl, lmask);
                                if hash_l <= hash_bound {
                                    seq_hashes.push(hash_l);
                                    //pos_to_seq_coord.push(i - l + 1);
                                    if !params.use_hpc {pos_to_seq_coord.push(tup.1[i-l+1]);} //if not HPCd need raw sequence positions
                                    else {pos_to_seq_coord.push(i-l+1);} //already HPCd so positions are the same
                                }
                            }
                        }
                    }
                }
                else { //kminmer
                    let yl : u64 = match xl[0] < xl[1] {
                        true => xl[0],
                        false => xl[1]
                    };
                    let hash_l = hash(yl, lmask);
                    if hash_l <= hash_bound {
                        seq_hashes.push(hash_l);

                        if !params.use_hpc {pos_to_seq_coord.push(tup.1[i-l+1]);} //if not HPCd need raw sequence positions
                        else {pos_to_seq_coord.push(i-l+1);} //already HPCd so positions are the same
                        //pos_to_seq_coord.push(i - l + 1);
                    }
                }
            } else {
                qs_min_val = u64::max_value();
                qs_min_pos = -1;
                lp = 0; xs = [0; 2]; xl = [0; 2];
                qs_size = 0;
                qs.clear();
                qs_pos.clear();
            }
        }
        let read_minimizers = Vec::<String>::new(); // unused everywhere it seems
        Read {id: inp_id.to_string(), minimizers: read_minimizers, minimizers_pos: pos_to_seq_coord, transformed: seq_hashes, seq: inp_seq_raw, corrected: false}
    } 

    pub fn label(&self, read_seq: String, read_minimizers: Vec<String>, read_minimizers_pos: Vec<usize>, read_transformed: Vec<u64>, corrected_map: &mut CorrMap) {
        corrected_map.insert(self.id.to_string(), (read_seq, read_minimizers, read_minimizers_pos, read_transformed));
    }

    pub fn read_to_kmers(&mut self, params: &Params) -> Repr {
        let k = params.k;
        let l = params.l;
        let mut output : Repr = Vec::new();
        for i in 0..(self.transformed.len()- k + 1) {
            let mut node : Kmer = Kmer::make_from(&self.transformed[i..i+k]);
            let mut seq_reversed = false;
            if REVCOMP_AWARE { 
                let (node_norm, reversed) = node.normalize(); 
                node = node_norm;
                seq_reversed = reversed;
            } 
            let seq = self.seq[self.minimizers_pos[i] as usize..(self.minimizers_pos[i+k-1] as usize + l)].to_string();
            //if seq_reversed {
            //    seq = utils::revcomp(&seq);
            //}

            /* // TODO incorporate that code somehow into writing an adequate sequence 
              // not just the first sequence that appears for a kmer (at abundance 2)
              // but rather the seq that has closest length to median length for that kmer
               let mut inserted = false;
               if kmer_seqs.contains_key(&node) { 
               let median_seq_len : usize = median(kmer_seqs_lens.get(&node).unwrap()) as usize;
            //println!("node: {} seqlen: {}",node.print_as_string(),seq.len());
            // insert that sequence if it's closer to the median length than the current
            // inserted string
            if ((seq.len() - median_seq_len) as f64).abs() < ((kmer_seqs.get(&node).unwrap().len() - median_seq_len) as f64).abs()
            { 
            kmer_seqs.insert(node.clone(), seq.clone());
            inserted = true;
            }
            }
            else
            {
            kmer_seqs.insert(node.clone(), seq.clone());
            inserted = true;
            }
            kmer_seqs_lens.entry(node.clone()).or_insert(Vec::new()).push(seq.len() as u32);
            */

            let origin = "*".to_string(); // uncomment the line below to track where the kmer is coming from (but not needed in production)
            //let origin = format!("{}_{}_{}", self.id, self.minimizers_pos[i].to_string(), self.minimizers_pos[i+k-1].to_string()); 

            let position_of_second_minimizer = match seq_reversed {
                true => self.minimizers_pos[i+k-1] - self.minimizers_pos[i+k-2],
                false => self.minimizers_pos[i+1] - self.minimizers_pos[i]
            };
            let position_of_second_to_last_minimizer = match seq_reversed {
                true => self.minimizers_pos[i+1] - self.minimizers_pos[i],
                false => self.minimizers_pos[i+k-1] - self.minimizers_pos[i+k-2]
            };
            let shift = (position_of_second_minimizer, position_of_second_to_last_minimizer);
            output.push((node, seq, seq_reversed, origin, shift));
        }
        output
    }
    pub fn poa_correct(&mut self, int_to_minimizer: &HashMap<u64, String>, poa_map: &mut HashMap<String, Vec<String>>, buckets: &Buckets, params : &Params, mut corrected_map: &mut CorrMap, reads_by_id: &HashMap<String, Read>) {

	// poa alignment scoring parameters
	let score = |a: u64, b: u64| if a == b {1i32} else {-1i32};
        let scoring = poa::Scoring::new(-1, -1, score);
        //let scoring = poa::Scoring::new(-1, 0, score); // the hope is that a 0 gap extend penalty somehow approximates semiglobal, but that's not quite true
        // other alignment parameters
        let dist_threshold = 0.15; // mash distance cut-off for read recruitment
        //let top_x_aligned_reads = 0; // get the 10 best read alignments per template
        let poa_global_min_score = std::i32::MIN; // discard all alignments below that score (discarded when top_X_aligned_read > 0)
        //let mut poa_global_min_score = -10; 
        let debug = params.debug;
        let n = params.n;
        let l = params.l;
        let read_minimizers_pos = &self.minimizers_pos;
        let read_transformed = &self.transformed;
        let seq_id = &self.id;
        let seq_str = &self.seq;
        let mut added_reads : HashSet<String> = HashSet::new();
        let mut bucket_reads = Vec::<&Read>::new();
        let mut poa_ids = Vec::<String>::new();
        let mut aligner = poa::Aligner::new(scoring, &read_transformed, Some(seq_str), Some(read_minimizers_pos));
        // populate bucket_reads with reads that share n consecutive minimizers with template
        added_reads.insert(self.id.clone());  
        for i in 0..read_transformed.len()-n+1 {
            let bucket_idx = &utils::normalize_vec(&read_transformed[i..i+n].to_vec());
            let entry = &buckets[bucket_idx];
            for id in entry.iter() {
                let query = &reads_by_id[id];
                if !added_reads.contains(&query.id) {
                    bucket_reads.push(query);
                    added_reads.insert(query.id.clone());  
                }
            }   
        }
        // filter bucket_reads so that it only contains reads below a mash distance of template
        let mut bucket_reads : Vec<(&Read, f64)> = bucket_reads.iter().map(|seq| (*seq, minimizers::dist(self, seq, &params))).filter(|(_seq, dist)| *dist < dist_threshold).collect();
        // sort reads by their distance
        bucket_reads.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let max_poa_reads = 80;
        if bucket_reads.len() > max_poa_reads {
            bucket_reads = bucket_reads[..max_poa_reads].to_vec(); // limit to 50 lowest-distance reads per POA correction
        }
        //if bucket_reads.len() == 0 {println!("Read has no neighbors");}
        // do a first pass aligning reads to the template to get only the top X best scoring reads
        // TODO this code isn't up to date as it only aligns in forward direction
        /*
        let mut alignment_scores = Vec::new();
        if top_X_aligned_reads > 0 
        {
            for i in 0..bucket_reads.len() {
                aligner.global(&bucket_reads[i].0.transformed);
                let score = aligner.alignment().score;
                alignment_scores.push((score,i));
            }
            // sort alignment scores decreasing
            alignment_scores.sort_by(|a, b| b.cmp(a));
            //if debug { println!("read alignment scores: {:?}",alignment_scores); }

            // threshold becomes the X-th highest alignment
            if alignment_scores.len() > top_X_aligned_reads
            {
                poa_global_min_score = alignment_scores[top_X_aligned_reads].0;
            }
            else{
                poa_global_min_score = std::i32::MIN;
            }
        }
        */
        let mut nb_aln_forward = 0;
        let mut nb_aln_backward = 0;
        for bucket_read in &bucket_reads {
            poa_ids.push(bucket_read.0.id.to_string());
            // don't know how to save alignments so i'll just waste time and recompute the best
            // scoring alignment.
            // it's a TODO optimization
            let read = &bucket_read.0;
            let seq = &read.seq;
            let pos = &read.minimizers_pos;
            aligner.semiglobal(&read.transformed, Some(seq), Some(pos));
            let fwd_score = aligner.alignment().score;
            if debug {println!("--- Forward alignment score: {} (ID: {})\nMinimizer-space: {}\n{}\n---", aligner.alignment().score, read.id.to_string(), pretty_minvec(&read.transformed), aligner.print_aln());}
            let mut rev_read = read.transformed.clone();
            rev_read.reverse();
            let rev_seq = utils::revcomp(&seq);
            let mut rev_minim_pos = pos.clone();
            rev_minim_pos.reverse();
            for pos in &mut rev_minim_pos {
                *pos = seq.len() - l - *pos;
            }
            aligner.semiglobal(&rev_read, Some(&rev_seq), Some(&rev_minim_pos));
            if debug { println!("--- Backward alignment score: {} (ID: {})\nMinimizer-space: {}\n{}\n---", aligner.alignment().score, read.id.to_string(), pretty_minvec(&rev_read), aligner.print_aln());}
            let bwd_score = aligner.alignment().score;
            //let mut aln_ori = "";
            if std::cmp::max(fwd_score, bwd_score) >= poa_global_min_score { 
               if fwd_score > bwd_score { 
                    aligner.semiglobal(&read.transformed, Some(seq), Some(pos));
                    //aln_ori = "fwd";
                    nb_aln_forward += 1;
                } else { 
                    aligner.semiglobal(&rev_read, Some(&rev_seq), Some(&rev_minim_pos));
                    //aln_ori = "bwd";
                    nb_aln_backward += 1;
                }
                aligner.add_to_graph(); 
            }
            //if debug { aligner.traceback.print(aligner.graph(), bucket_reads[i].0.transformed.clone()); } // prints a confusing traceback (and crashes)
        }
        let (consensus, consensus_edge_seqs) = aligner.poa.consensus(&params);
        //println!("consensus()/consensus_edge_seqs() lens: {} / {}", consensus.len(), consensus_edge_seqs.len());
        let (consensus, consensus_edge_seqs) = aligner.consensus_boundary(&consensus, &consensus_edge_seqs, &read_transformed, debug);
        let consensus_read = consensus.iter().map(|minim| int_to_minimizer[minim].to_string()).collect::<Vec<String>>();
        let len_before_poa  = self.transformed.len();
        let len_after_poa  = consensus_read.len();
        if debug { println!("Length of template before/after POA: {} / {} (ID: {})", len_before_poa, len_after_poa, seq_id);}
        if debug { println!("Number of bucketed reads aligned forwards/backwards: {} / {} ", nb_aln_forward, nb_aln_backward);}
        let mut consensus_str = String::new();
        let mut pos_idx = 0;
        let mut consensus_pos = Vec::<usize>::new();
        if consensus.is_empty() {return}
        for insert in consensus_edge_seqs.iter() {
            consensus_pos.push(pos_idx);
            consensus_str.push_str(insert);
            pos_idx += insert.len();
        }
        consensus_pos.push(pos_idx);
        consensus_str.push_str(&int_to_minimizer[&consensus[consensus.len()-1]]);
        let mut corrected_count = 0;
        let mut threshold = params.correction_threshold;
        if params.correction_threshold == 0 {threshold = 0 as i32;}
        for (read, _dist) in bucket_reads.iter_mut() {
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
