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
use fasthash::city;
use itertools::Itertools;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::iter::FromIterator;
use xx_bloomfilter::Bloom;
use super::RacyBloom;
use super::read::Read;
use nthash::{ntc64, NtHashIterator};
const lmer_frequency_based : bool = false;

// Computes the distance between two reads which are already in minimizer-space.

pub fn dist(temp: &Read, other: &Read, params: &Params) -> f64 {
    let s1_set: HashSet<_> = HashSet::from_iter(temp.transformed.iter());
    let s2_set: HashSet<_> = HashSet::from_iter(other.transformed.iter());
    let inter: HashSet<_> = s1_set.intersection(&s2_set).collect();
    let union: HashSet<_> = s1_set.union(&s2_set).collect();
    let distance = params.distance;
    match distance {
        0 => {return 1.0 - ((inter.len() as f64) / (union.len() as f64))}
        1 => {return 1.0 - ((inter.len() as f64) / (s1_set.len() as f64))}
        2 => {
            let jaccard = (inter.len() as f64) / (union.len() as f64);
            let mash : f64 = -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / params.l as f64;
            return mash
        }
        _ => {
            let jaccard = (inter.len() as f64) / (union.len() as f64);
            let mash : f64 = (-1.0 / params.k as f64) * ((2.0 * jaccard) / (1.0 + jaccard)).ln();
            return mash
        }
    }
}
pub fn normalize_minimizer(lmer: &str) -> String {
    let mut res = lmer.to_string();
    if revcomp_aware {
        let rev = utils::revcomp(&lmer);
        if rev < res {res = rev.clone();}
    }
    res
}

// https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter
pub fn minimizers_preparation(mut params: &mut Params, filename: &PathBuf, file_size: u64, lmer_counts: &HashMap<String, u32>) -> (HashMap<String, u64>, HashMap<u64, String>) {
    let l = params.l;
    let density = params.density;
    let mut list_minimizers : Vec<String> = Vec::new();
    let mut count_vec : Vec<(&String, &u32)> = lmer_counts.into_iter().collect();
    let mut max_threshold = params.lmer_counts_max;
    let mut min_threshold = params.lmer_counts_min;
    let mut skip : HashSet<String> = HashSet::new();
    // the following code replaces what i had before:
    // https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter
    let multi_prod = (0..l).map(|i| vec!('A', 'C', 'T', 'G')).multi_cartesian_product();
    for lmer_vec in multi_prod {
        let lmer : String = lmer_vec.into_iter().collect();
        if revcomp_aware {
            let lmer_rev = utils::revcomp(&lmer);
            if lmer > lmer_rev {continue;} // skip if not canonical
        }
        list_minimizers.push(lmer.to_string());
    }
    count_vec.iter().filter(|tup| tup.1 >= &max_threshold || tup.1 <= &min_threshold).map(|tup| tup.0.to_string()).collect::<Vec<String>>().iter().for_each(|x| {skip.insert(x.to_string());});
    
    let mut minimizer_to_int : HashMap<String, u64> = HashMap::new();
    let mut int_to_minimizer : HashMap<u64, String> = HashMap::new();
    let mut skips = 0;
        // assign numbers to minimizers, the regular way
        for lmer in list_minimizers {
            let mut hash = (ntc64(lmer.as_bytes(), 0, l)) as u64;
            let mut hash_new = hash as f64;
            hash_new = (hash_new) / (u64::max_value() as f64);
            if skip.contains(&lmer) { 
                hash_new = hash_new.sqrt().sqrt().sqrt();
                skips += 1;
            }
            if (hash_new as f64) <= (density as f64) {
                minimizer_to_int.insert(lmer.to_string(), hash);
                int_to_minimizer.insert(hash, lmer.to_string());
                // also insert the info for the revcomp minimizer, will allow to avoid normalizing
                // later 
                let lmer_rev = utils::revcomp(&lmer);
                minimizer_to_int.insert(lmer_rev.to_string(), hash);
                int_to_minimizer.insert(hash, lmer_rev.to_string());
            }
        }
    
    println!("Selected {} minimizer ID's and {} sequences.",int_to_minimizer.len(), minimizer_to_int.len());
    println!("{} frequent l-mers skipped.", skips);
    (minimizer_to_int, int_to_minimizer)
}

pub fn uhs_preparation(mut params: &mut Params, uhs_filename: &str) -> RacyBloom {
    let l = params.l;
    let mut uhs_bloom : RacyBloom = RacyBloom::new(Bloom::new(if params.use_bf {500_000_000} else {1}, 1_000_000_000_000_000));
    let ufile = match File::open(uhs_filename) {
        Ok(ufile) => ufile,
        Err(error) => panic!("Error opening UHS file: {:?}.", error),
    };
    let mut ureader = BufReader::new(ufile);
    let mut umer = String::new();
    loop {
        match ureader.read_line(&mut umer) {
            Err(_) | Ok(0) => break,
            Ok(_) => {
                let hash = ntc64(umer.as_bytes(), 0, l);
                uhs_bloom.get().check_and_add(&hash);
                umer.clear();
            }
        }
    }
    println!("All universal k-mers read.");
    uhs_bloom
}

pub fn lcp_preparation(mut params: &mut Params, lcp_filename: &str) -> RacyBloom {
    let l = params.l;
    let mut lcp_bloom : RacyBloom = RacyBloom::new(Bloom::new(if params.use_bf {500_000_000} else {1}, 1_000_000_000_000_000));
    let lcpfile = match File::open(lcp_filename) {
        Ok(lcpfile) => lcpfile,
        Err(error) => panic!("Error opening LCP file: {:?}.", error),
    };
    let mut lcp_reader = BufReader::new(lcpfile);
    let mut core = String::new();
    loop {
        match lcp_reader.read_line(&mut core) {
            Err(_) | Ok(0) => break,
            Ok(_) => {
                if core.len() == l {
                    let hash = ntc64(core.as_bytes(), 0, l);
                    lcp_bloom.get().check_and_add(&core);
                    core.clear();
                }
            }
        }
    }
    println!("All LCP core substrings read.");
    lcp_bloom
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
