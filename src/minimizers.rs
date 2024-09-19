//use std::i64;
use super::utils;
use super::Params;
use super::REVCOMP_AWARE;
use std::collections::HashMap;
use std::collections::HashSet;
//use pbr::ProgressBar;
//use bio::io::fasta;
//use std::path::PathBuf;
use itertools::Itertools;
use std::fs::File;
use std::io::{BufRead, BufReader};
//use std::path::Path;
use std::iter::FromIterator;
use xx_bloomfilter::Bloom;
use super::RacyBloom;
use super::read::Read;
use nthash::{ntc64};//, NtHashIterator};

// Computes the distance between two reads which are already in minimizer-space.

pub fn dist(temp: &Read, other: &Read, params: &Params) -> f64 {
    let s1_set: HashSet<_> = HashSet::from_iter(temp.transformed.iter());
    let s2_set: HashSet<_> = HashSet::from_iter(other.transformed.iter());
    let inter: HashSet<_> = s1_set.intersection(&s2_set).collect();
    let union: HashSet<_> = s1_set.union(&s2_set).collect();
    let distance = params.distance;
    match distance {
        0 => {1.0 - ((inter.len() as f64) / (union.len() as f64))}
        1 => {1.0 - ((inter.len() as f64) / (s1_set.len() as f64))}
        2 => {
            let jaccard = (inter.len() as f64) / (union.len() as f64);
            let mash : f64 = -1.0 * ((2.0 * jaccard) / (1.0 + jaccard)).ln() / params.l as f64;
            mash
        }
        _ => {
            let jaccard = (inter.len() as f64) / (union.len() as f64);
            let mash : f64 = (-1.0 / params.k as f64) * ((2.0 * jaccard) / (1.0 + jaccard)).ln();
            mash
        }
    }
}
/*pub fn normalize_minimizer(lmer: &str) -> String {
    let mut res = lmer.to_string();
    if revcomp_aware {
        let rev = utils::revcomp(&lmer);
        if rev < res {res = rev.clone();}
    }
    res
}*/

// https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter
pub fn minimizers_preparation(params: &mut Params, lmer_counts: &HashMap<String, u32>) -> (HashMap<String, u64>, HashMap<u64, String>) {
    let l = params.l;
    let density = params.density;
    let mut list_minimizers : Vec<String> = Vec::new();
    let count_vec : Vec<(&String, &u32)> = lmer_counts.iter().collect();
    let max_threshold = params.lmer_counts_max;
    let min_threshold = params.lmer_counts_min;
    let mut skip : HashSet<String> = HashSet::new();
    // large l-mer fix: use lmer counts to enumerate them, instead of cartesian product
    if lmer_counts.len() > 0 {
        count_vec.iter().map(|tup| tup.0).collect::<Vec<&String>>().iter().for_each(|x| {
            list_minimizers.push(std::cmp::min(x.to_string(),utils::revcomp(x)))});
    }
    else
    {
        // the following code replaces what i had before:
        // https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter
        let multi_prod = (0..l).map(|_i| vec!('A', 'C', 'T', 'G')).multi_cartesian_product();
        for lmer_vec in multi_prod {
            let lmer : String = lmer_vec.into_iter().collect();
            if REVCOMP_AWARE {
                let lmer_rev = utils::revcomp(&lmer);
                if lmer > lmer_rev {continue;} // skip if not canonical
            }
            list_minimizers.push(lmer.to_string());
        }
    }
    count_vec.iter().filter(|tup| tup.1 >= &max_threshold || tup.1 <= &min_threshold).map(|tup| tup.0.to_string()).collect::<Vec<String>>().iter().for_each(|x| {
        skip.insert(x.to_string());
        skip.insert(utils::revcomp(x));
    });
    
    let mut minimizer_to_int : HashMap<String, u64> = HashMap::new();
    let mut int_to_minimizer : HashMap<u64, String> = HashMap::new();
    let mut skips = 0;
    // assign numbers to minimizers, the regular way
    for lmer in list_minimizers {
        let hash = (ntc64(lmer.as_bytes(), 0, l)) as u64;
        let mut hash_new = hash as f64;
        hash_new /= u64::max_value() as f64;
        if skip.contains(&lmer) { 
            hash_new = 1.0;
            //println!("{} skipping lmer", lmer);
            skips += 1;
        }
        if (hash_new as f64) <= (density as f64) {
            //println!("{} lmer {} hash {} density", lmer,hash_new,density);
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

pub fn uhs_preparation(params: &mut Params, uhs_filename: &str) -> RacyBloom {
    let l = params.l;
    let uhs_bloom : RacyBloom = RacyBloom::new(Bloom::new(if params.use_bf {500_000_000} else {1}, 1_000_000_000_000_000));
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

pub fn lcp_preparation(params: &mut Params, lcp_filename: &str) -> RacyBloom {
    let l = params.l;
    let lcp_bloom : RacyBloom = RacyBloom::new(Bloom::new(if params.use_bf {500_000_000} else {1}, 1_000_000_000_000_000));
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
                    //let hash = ntc64(core.as_bytes(), 0, l);
                    lcp_bloom.get().check_and_add(&core);
                    core.clear();
                }
            }
        }
    }
    println!("All LCP core substrings read.");
    lcp_bloom
}
