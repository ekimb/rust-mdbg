//mod super::params;
use std::i64;
use super::utils;
use super::Params;
use rayon::prelude::*;
use super::revcomp_aware;
use std::collections::HashMap;
use std::collections::HashSet;
use pbr::ProgressBar;
use bio::io::fasta;
use std::path::PathBuf;
use strsim::levenshtein;
use fasthash::city;
use itertools::Itertools;

use nthash::{ntc64, NtHashIterator};

const lmer_frequency_based : bool = false;


/* why we need levenshtein balls as minimizers?

there are 13147 kmers in the reference genome graph of dmel chr4, k10-p0.01-l12

in reads with error rate 0.0075:
43114 kmers in graph-k10-p0.01-l12.sequences
number of kmers from genomegraph-k10-p0.01-l12.sequences that are in graph-k10-p0.01-l12.sequences 13112 (99.73)% , 35 are not

in reads with erro rate 0.02:

# python ../rust-mhdbg/utils/compare_kmers.py genomegraph-0.02-k10-p0.01.sequences graph-0.02-k10-p0.01.sequences:w
15244 kmers in graph-0.02-k10-p0.01.sequences
number of kmers from genomegraph-k10-p0.01-l12.sequences that are in graph-0.02-k10-p0.01.sequences 8951 (68.08)% , 4196 are not

bottom line: with k=10, couldn't retrieve enough solid (minabund>=2) kmers (at read error rate is 2%)
*/

pub fn normalize_minimizer(lmer: &String) -> String
{
    let mut res = lmer.clone();
    if revcomp_aware {
        let rev = utils::revcomp(&lmer);
        if rev < res { res = rev; }
    }
    res
}

pub fn minhash(seq: String, params: &Params, int_to_minimizer : &HashMap<u64, String>) -> (Vec<String>, Vec<u32>, Vec<u64>)
{
    let size_miniverse = params.size_miniverse as u64;
    let density = params.density;
    let l = params.l;
    let mut read_minimizers = Vec::<String>::new();
    let mut read_minimizers_pos = Vec::<u32>::new();
    let mut read_transformed = Vec::<u64>::new();
    for i in 0..seq.len()-l+1 {
        let lmer = &seq[i..i+l];
        let hash = ntc64(lmer.as_bytes(), 0, l);
        if (hash != 0) && (hash as f64) < u64::max_value() as f64 * density/(l as f64) {
            read_minimizers.push(normalize_minimizer(&lmer.to_string()));
            read_minimizers_pos.push(i as u32);
            read_transformed.push(hash);
        }
    }
    
    //let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
    // convert minimizers to their integer representation

    (read_minimizers, read_minimizers_pos, read_transformed)
}


// https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter

pub fn minimizers_preparation(mut params: &mut Params, filename :&PathBuf, file_size: u64, levenshtein_minimizers: usize) -> (HashMap<String,u64>, HashMap<u64,String>) {

    let l = params.l;
    let mut list_minimizers : Vec<String> = Vec::new();
    
    // the following code replaces what i had before:
    // https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter
    let multi_prod = (0..l).map(|i| vec!('A','C','T','G'))
            .multi_cartesian_product();

//    for lmer in kproduct("ACTG".to_string(), l as u32) {
    for lmer_vec in multi_prod {
        let lmer :String = lmer_vec.into_iter().collect();
        //println!("testing minimizer {}",lmer.to_string());
        //println!("found minimizer {}",lmer.to_string());
        if revcomp_aware {
            let lmer_rev = utils::revcomp(&lmer);
            if lmer > lmer_rev {continue;} // skip if not canonical
       }
        list_minimizers.push(lmer);
    }
   
    let mut minimizer_to_int : HashMap<String,u64> = HashMap::new();
    let mut int_to_minimizer : HashMap<u64,String> = HashMap::new();
    let mut minim_idx : u32 = 0;
        // assign numbers to minimizers, the regular way
        for lmer in list_minimizers
        {
            let hash = (ntc64(lmer.as_bytes(), 0, l)) as u64;
            minimizer_to_int.insert(lmer.to_string(),  hash);
            int_to_minimizer.insert(hash,         lmer.to_string());
            minim_idx += 1;
        }
    
    println!("selected {} minimizer ID's, {} sequences",int_to_minimizer.len(), minimizer_to_int.len());
    (minimizer_to_int, int_to_minimizer)
}

