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
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::iter::FromIterator;
use super::read::Read;
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
pub fn dist(temp: &Read, other: &Read, params: &Params) -> f64 {
    let s1_set: HashSet<_> = HashSet::from_iter(temp.transformed.iter());
    let s2_set: HashSet<_> = HashSet::from_iter(other.transformed.iter());
    let inter: HashSet<_> = s1_set.intersection(&s2_set).collect();
    let union: HashSet<_> = s1_set.union(&s2_set).collect();
    let distance = params.distance;
    match distance {
        0 => {
            return 1.0 - ((inter.len() as f64) / (union.len() as f64))
        }
        1 => {
            return 1.0 - ((inter.len() as f64) / (s1_set.len() as f64))
        }
        2 => {
            let jaccard = (inter.len() as f64) / (union.len() as f64);
            let mash: f64 = -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / params.l as f64;
            return mash
        }
        _ => {
            let jaccard = (inter.len() as f64) / (union.len() as f64);
            let mash: f64 = (-1.0 / params.k as f64) * ((2.0 * jaccard) / (1.0 + jaccard)).ln();
            return mash
        }
    }
}
pub fn lmer_counting(lmer_counts: &mut HashMap<String,u32>, filename :&PathBuf, file_size :u64, params: &mut Params) {
    let l = params.l;
    let mut pb = ProgressBar::new(file_size);
    let reader = fasta::Reader::from_file(&filename).unwrap();

    for record in reader.records() {
        let record = record.unwrap();
        let seq = record.seq();

        let seq_str = String::from_utf8_lossy(seq);
        pb.add(record.seq().len() as u64 + record.id().len() as u64); // approx size of entry

        if seq_str.len() < l { continue; }
    
        //let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
        for start in 0..(&seq_str.len()-l+1)
        {
            let mut lmer = String::from(&seq_str[start..start+l]);
            if revcomp_aware {
                let lmer_rev = utils::revcomp(&lmer);
                lmer = std::cmp::min(lmer, lmer_rev);
            }
            let count = lmer_counts.entry(lmer).or_insert(0);
            *count += 1;
        }
    }
    params.average_lmer_count = lmer_counts.values().sum::<u32>() as f64 / lmer_counts.len() as f64;
    println!("average lmer count: {}",params.average_lmer_count);
}
pub fn normalize_minimizer(lmer: &String) -> String
{
    let mut res = lmer.clone();
    if revcomp_aware {
        let rev = utils::revcomp(&lmer);
        if rev < res { res = rev; }
    }
    res
}

pub fn minhash(seq: String, params: &Params, int_to_minimizer : &HashMap<u64, String>, lmer_counts: &mut HashMap<String, u32>) -> (Vec<String>, Vec<u32>, Vec<u64>)
{
    let size_miniverse = params.size_miniverse as u64;
    let density = params.density;
    let l = params.l;
    let mut read_minimizers = Vec::<String>::new();
    let mut read_minimizers_pos = Vec::<u32>::new();
    let mut read_transformed = Vec::<u64>::new();
    for i in 0..seq.len()-l+1 {
        let mut lmer = &seq[i..i+l];
        let lmer_rev = utils::revcomp(&lmer);
        lmer = std::cmp::min(lmer, &lmer_rev);
        let hash = ntc64(lmer.as_bytes(), 0, l);
        if (hash != 0) && (hash as f64) < u64::max_value() as f64 * density/(l as f64) {
            read_minimizers.push(lmer.to_string());
            read_minimizers_pos.push(i as u32);
            read_transformed.push(hash);
        }
    }
    
    //let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
    // convert minimizers to their integer representation

    (read_minimizers, read_minimizers_pos, read_transformed)
}

pub fn minhash_uhs(seq: String, params: &Params, int_to_minimizer : &HashMap<u64, String>, lmer_counts: &mut HashMap<String, u32>, uhs_kmers: &HashMap<String, u32>) -> (Vec<String>, Vec<u32>, Vec<u64>)
{
    let size_miniverse = params.size_miniverse as u64;
    let density = params.density;
    let l = params.l;
    let mut read_minimizers = Vec::<String>::new();
    let mut read_minimizers_pos = Vec::<u32>::new();
    let mut read_transformed = Vec::<u64>::new();
    for i in 0..seq.len()-l+1 {
        let mut lmer = &seq[i..i+l];
        let mut lmer = normalize_minimizer(&lmer.to_string());
        let hash = ntc64(lmer.as_bytes(), 0, l);
        if (hash != 0) && (hash as f64) < u64::max_value() as f64 * density/(l as f64) && uhs_kmers[&lmer] == 1 {
            read_minimizers.push(lmer.to_string());
            read_minimizers_pos.push(i as u32);
            read_transformed.push(hash);
        }
    }
    
    //let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
    // convert minimizers to their integer representation

    (read_minimizers, read_minimizers_pos, read_transformed)
}
pub fn minhash_window(seq: String, params: &Params, int_to_minimizer : &HashMap<u64, String>, minimizer_to_int : &HashMap<String, u64>, lmer_counts: &mut HashMap<String, u32>) -> (Vec<String>, Vec<u32>, Vec<u64>)
{
    let l = params.l;
    let w = params.w;
    let mut read_minimizers = Vec::<String>::new();
    let mut read_minimizers_pos = Vec::<u32>::new();
    let mut read_transformed = Vec::<u64>::new();
    for i in 0..seq.len()-w+1 {
        let window = &seq[i..i+w];
        let mut min = String::new();
        let mut min_hash = u64::max_value();
        let mut min_pos = usize::max_value();
        for j in 0..w-l+1 {
            let mut lmer = &window[j..j+l];
            let mut lmer = normalize_minimizer(&lmer.to_string());
            let hash = minimizer_to_int[&lmer];
            if !lmer.contains("N") && hash < min_hash {
                min = lmer.to_string();
                min_hash = hash;
                min_pos = i+j;
            }
        }
        if (read_minimizers_pos.len() != 0) && min_pos != usize::max_value() {
            if min_pos as u32 != read_minimizers_pos[read_minimizers_pos.len()-1] {
                read_minimizers.push(min.to_string());
                read_minimizers_pos.push(min_pos as u32);
                read_transformed.push(min_hash);
            }
        }
        else if min_pos != usize::max_value() {
            read_minimizers.push(min.to_string());
            read_minimizers_pos.push(min_pos as u32);
            read_transformed.push(min_hash);
        }
    }
    //let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
    // convert minimizers to their integer representation

    (read_minimizers, read_minimizers_pos, read_transformed)
}


// https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter

pub fn minimizers_preparation(mut params: &mut Params, filename :&PathBuf, file_size: u64, levenshtein_minimizers: usize, lmer_counts: &HashMap<String, u32>) -> (HashMap<String,u64>, HashMap<u64,String>) {

    let l = params.l;
    let density = params.density;
    let mut list_minimizers : Vec<String> = Vec::new();
    let mut count_vec: Vec<(&String, &u32)> = lmer_counts.into_iter().collect();
    let mut threshold = 100;
    let mut skip : HashMap<String, bool> = HashMap::new();
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
        list_minimizers.push(lmer.to_string());
        skip.insert(lmer, false);
    }
    count_vec.iter().filter(|tup| tup.1 >= &threshold).map(|tup| tup.0.to_string()).collect::<Vec<String>>().iter().for_each(|x| *skip.entry(x.to_string()).or_insert(false) = true);
    
    let mut minimizer_to_int : HashMap<String,u64> = HashMap::new();
    let mut int_to_minimizer : HashMap<u64,String> = HashMap::new();
    let mut minim_idx : u32 = 0;
    let mut skips = 0;
        // assign numbers to minimizers, the regular way
        for lmer in list_minimizers
        {
            let mut hash = (ntc64(lmer.as_bytes(), 0, l)) as u64;
            if skip[&lmer] {
                hash = u64::max_value();
                skips += 1;
            }
            minimizer_to_int.insert(lmer.to_string(),  hash);
            int_to_minimizer.insert(hash,         lmer.to_string());
            minim_idx += 1;
        }
    
    println!("selected {} minimizer ID's, {} sequences",int_to_minimizer.len(), minimizer_to_int.len());
    println!("{} frequent l-mers skipped", skips);
    (minimizer_to_int, int_to_minimizer)
}

pub fn uhs_preparation(mut params: &mut Params, uhs_filename : &str) -> HashMap<String, u32> {

    let l = params.l;
    let mut uhs_kmers = HashMap::<String, u32>::new();
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
        uhs_kmers.insert(lmer, 0);
    }
    let mut count = 0;
    if let Ok(lines) = read_lines(uhs_filename) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines {
            if let Ok(mut ip) = line {
                *uhs_kmers.entry(normalize_minimizer(&ip)).or_insert(0) = 1;
                count += 1;
            }
        }
    }
    println!("selected {} umers", count);
    uhs_kmers
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}