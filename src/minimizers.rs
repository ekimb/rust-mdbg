//mod super::params;
use super::utils;
use super::Params;
use super::revcomp_aware;
use std::collections::HashMap;
use std::collections::HashSet;
use pbr::ProgressBar;
use bio::io::fasta;
use std::path::PathBuf;
use strsim::levenshtein;
use fasthash::city;
use itertools::Itertools;
use std::borrow::Cow;


const lmer_frequency_based : bool = true;
const simple_correct : bool = false;
pub fn correct_kmers(lmer_counts: &mut HashMap<String,u32>, seq_str : &str, params: &Params) -> String{
    let l = params.l;
    let threshold = params.t;
    let mut new_seq = seq_str.to_string();
    let mut errSwitch : bool = false;
    let mut errStartPos = 0;
    let mut errSuffPos = 0;
    let mut errEndPos = 0;
    if simple_correct {
        
        for start in 0..(&new_seq.len()-l+1) {
            let mut lmer = String::from(&new_seq[start..start+l]);
                if revcomp_aware {
                    let lmer_rev = utils::revcomp(&lmer);
                    lmer = std::cmp::min(lmer, lmer_rev);
                }
            //println!("{}", &seq_str[start..start+k]);
            let mut count = 0;
            //println!("{}", lmer);
            if lmer_counts.contains_key(&lmer) {
                count = lmer_counts[&lmer];
                //println!("K-mer count: {}", &count);
            }
            if count <= threshold {
                let levBall = levenshtein_ball(&lmer, 1);
                let mut maxNeighbor = String::new();
                let mut neighbor = String::new();
                let mut neighborMin = String::new();
                let mut maxCount = 0;

                for neighbor in levBall {
                    if revcomp_aware {
                        let neighborRev = utils::revcomp(&neighbor);
                        neighborMin = std::cmp::min(neighbor, neighborRev);
                    }
                    //println!("{}", old);
                    //println!("{}", neighbor);
                    //assert_eq!(neighbor.len(), l);
                    //println!("here");
                    if lmer_counts.contains_key(&neighborMin)  {
                        let mut neighborCount = lmer_counts[&neighborMin];
                        if neighborMin.len() > l {
                            neighborCount = lmer_counts[&neighborMin[0..l]];
                        }
                        if neighborCount > params.t {
                            maxCount = neighborCount;
                            maxNeighbor = neighborMin;
                            errSwitch = false;
                            //errStartPos = 0;
                            //errSuffPos = 0;
                            //errEndPos = 0;
                            break;
                        }
                    }

                }

                if maxCount <= threshold {
                    continue;
                }
                else {
                    new_seq = new_seq.replace(&lmer, &maxNeighbor);
                    

                }
            }
        }
        new_seq
    }
    else {
        let l = params.l;
        let threshold = params.t;
        let mut new_seq = seq_str.to_string();
        let mut errSwitch : bool = false;
        let mut errStartPos = 0;
        let mut errSuffPos = 0;
        let mut errEndPos = 0;
        for start in 0..(&new_seq.len()-l+1) {
            let mut lmer = String::from(&new_seq[start..start+l]);
                if revcomp_aware {
                    let lmer_rev = utils::revcomp(&lmer);
                    lmer = std::cmp::min(lmer, lmer_rev);
                }
            //println!("{}", &seq_str[start..start+k]);
            let mut count = 0;
            //println!("{}", lmer);
            if lmer_counts.contains_key(&lmer) {
                count = lmer_counts[&lmer];
                //println!("K-mer count: {}", &count);
            }
            else {
                //println!("K-mer count: {}", &count);
            }
            //println!("{}", count);
            if (count <= threshold && errSwitch) && (errSuffPos < errStartPos+l){
                //println!("Error continues");
                errSuffPos += 1;
                
            }
            else if (count <= threshold) && (!errSwitch) {
                //println!("Error starting");
                errStartPos = start;
                errSuffPos = start;
                errEndPos = start+l;
                errSwitch = true;
            }
            else if (count > threshold && errSwitch) {
                //println!("Error over");
                errSuffPos += 1;
                errSwitch = false;
                let mut old = String::from(&new_seq[start..start+l]);
                if revcomp_aware {
                    let oldRev = utils::revcomp(&old);
                    old = std::cmp::min(old, oldRev);
                }
            //println!("{}", &seq_str[start..start+k]);
                let mut count = 0;
                //println!("{}", old);
                if lmer_counts.contains_key(&old) {
                    count = lmer_counts[&lmer];
                    //println!("K-mer count: {}", &count);
                }
                if errSuffPos >= errEndPos {
                    let levBall = levenshtein_ball(&old, 1);
                    let mut maxCount = 0;
                    let mut maxNeighbor = String::new();
                    let mut neighbor = String::new();
                    let mut neighborMin = String::new();
                    let mut maxCount = 0;
                    let mut maxNeighbor = String::new();

                    for levNeighbor in levBall {
                        if revcomp_aware {
                            let neighborRev = utils::revcomp(&levNeighbor);
                            neighborMin = std::cmp::min(levNeighbor, neighborRev);
                        }
                        //println!("{}", old);
                        //println!("{}", neighbor);
                        //assert_eq!(neighbor.len(), l);
                        //println!("here");
                        if lmer_counts.contains_key(&neighborMin)  {
                            let mut neighborCount = lmer_counts[&neighborMin];
                            if neighborMin.len() > l {
                                neighborCount = lmer_counts[&neighborMin[0..l]];
                            }
                            if neighborCount > params.t {
                                neighborCount = maxCount;
                                maxNeighbor = neighborMin;
                                errSwitch = false;
                                //errStartPos = 0;
                                //errSuffPos = 0;
                                //errEndPos = 0;
                                break;
                            }
                        }

                    }

                    if maxCount <= threshold {
                        continue;
                    }
                    else {
                        new_seq = new_seq.replace(&old, &maxNeighbor);
                        let count = lmer_counts.entry(maxNeighbor).or_insert(0);
                        *count += 1;

                    }
                    //println!("Old l-mer: {}, new l-mer: {}", &old, &neighbor);
                    //errStartPos = 0;
                    //errSuffPos = 0;
                    //errEndPos = 0;
                    errSwitch = false;
            
                }
                else if errSuffPos != errEndPos {
                    let mut errSuffix = String::new();
                    errSuffix = String::from(&new_seq[errSuffPos..errEndPos]);
                    //println!("Suffix: {}", errSuffix);
                    //let suffixLen = errSuffix.len();
                    let levSuffBall = levenshtein_ball(&errSuffix, 1);
                    let mut maxCount = 0;
                    let mut maxNeighbor = String::new();
                    let mut neighborSuff = String::new();
                    for suff in levSuffBall {
                        //println!("LevSuff: {}", &suff);
                        //assert_eq!(neighbor.len(), l);
                        neighborSuff = String::from(&new_seq[errStartPos..errSuffPos]) + &suff;
                        if revcomp_aware {
                            let neighborRev = utils::revcomp(&neighborSuff);
                            neighborSuff = std::cmp::min(neighborSuff, neighborRev);
                        }
                        if lmer_counts.contains_key(&neighborSuff)  {
                            let neighborCount = lmer_counts[&neighborSuff];
                            //println!("Old count: {}", &count);
                            if neighborCount > params.t {
                                maxCount = neighborCount;
                                maxNeighbor = neighborSuff;
                                //new_seq = new_seq.replace(&old, &neighborSuff);
                                //println!("Old l-mer: {}, new l-mer: {}", &old, &neighborSuff);
                                errSwitch = false;
                                //errStartPos = 0;
                                //errSuffPos = 0;
                                //errEndPos = 0;
                                break;
                            }
                        }

                    }
                    if maxCount <= threshold {
                        continue;
                    }
                    else {
                        new_seq = new_seq.replace(&old, &maxNeighbor);
                        let count = lmer_counts.entry(maxNeighbor).or_insert(0);
                        *count += 1;

                    }
                    errSwitch = false;
                    //errStartPos = 0;
                    //errSuffPos = 0;
                    //errEndPos = 0;
                    errSuffix = String::new();
                    
                }
            
                //println!("maxcount: {}", maxCount);
                //if maxCount <= threshold {
                    //println!("Threshold not met");
                    //let levBallRev = levenshtein_ball(&kmer_rev, 1);
                    //let mut maxCount = 0;
                    //let mut maxNeighbor = String::new();
                    //for neighbor in levBallRev {
                        //if kmer_counts.contains_key(&neighbor)  {
                        // let neighborCount = kmer_counts[&neighbor];
                            //if neighborCount > maxCount {
                            // maxCount = neighborCount;
                            // maxNeighbor = utils::revcomp(&neighbor);
                        // }
                        //}
                // }
                    //if maxCount <= threshold {
                // continue;
                // }
                    
            // }
                //println!("Switched kmer with neighbor with count {}", maxCount);
                //let mut old = String::from(&seq_str[errStartPos..errEndPos]);
                
                //println!("{}", &seq_str[start..start+k]);

            }
        }
        new_seq
    }
    
}
pub fn print_hist(lmer_counts: &HashMap<String,u32>) {
    let mut hist = vec![0; 256];
    for (_lmer, count) in lmer_counts {
        let idx = if *count > 255 {
            255 as usize
        } else {
            *count as usize
        };
        hist[idx] += 1;
    }
    for i in 0..256 {
        println!("{}\t{}", i, hist[i]);
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

/* why we need levenshtein balls as minimizers?

there are 13147 kmers in the reference genome graph of dmel chr4, k10-p0.01-l12

in reads with error rate 0.0075:
43114 kmers in graph-k10-p0.01-l12.sequences
number of kmers from genomegraph-k10-p0.01-l12.sequences that are in graph-k10-p0.01-l12.sequences 13112 (99.73)% , 35 are not

in reads with erro rate 0.02:

# python ../rust-mhdbg/utils/compare_kmers.py genomegraph-0.02-k10-p0.01.sequences graph-0.02-k10-p0.01.sequences:w
15244 kmers in graph-0.02-k10-p0.01.sequences
number of kmers from genomegraph-k10-p0.01-l12.sequences that are in graph-0.02-k10-p0.01.sequences 8951 (68.08)% , 4196 are not

bottom line: with k=10, can't enough solid (minabund>=2) kmers when read error rate is 2%
*/

fn levenshtein_ball(lmer: &String, ball_size: usize) -> Vec<String>
{
    match ball_size
    {
        1 => levenshtein_ball_1(lmer),
        2 => levenshtein_ball_2(lmer),
        _ => panic!("unsupported levenshtein ball size"),
    }        
}

fn levenshtein_ball_1(lmer: &String) -> Vec<String>
{
    let mut res = HashSet::new();
    // all mutations
    for pos in 0..lmer.len()
    {
        for c in vec!('A','C','T','G')
        {
            if lmer.chars().nth(pos).unwrap() == c {continue;}
            let mut lmer_mutated = String::new();
            if pos > 0 {
                lmer_mutated.push_str(&lmer[..pos].to_string());
            }
            lmer_mutated.push_str(&c.to_string());
            if pos <= lmer.len()-1 {
                lmer_mutated.push_str(&lmer[pos+1..].to_string());
            }

            assert_eq!(levenshtein(&lmer_mutated,&lmer), 1);
            res.insert(lmer_mutated.to_string());
        }
    }
    // all deletions (yes, now, my minimizers will have length l and some l-1..)
    // except that we disallow deleting the first and last chars
    // otherwise we end up getting the adjacent minimizer too
    for pos in 1..(lmer.len()-1)
    {
        let mut lmer_mutated = lmer.clone();
        lmer_mutated.remove(pos);
        //println!("b4 {} pos {} a8 {}",&lmer,pos, &lmer_mutated);

        assert_eq!(levenshtein(&lmer_mutated,&lmer), 1);
        res.insert(lmer_mutated.to_string());
    }
    // all insertions
    // (now minimizer has length l+1 too..)
    for pos in 1..(lmer.len()-1)
    {
        for c in vec!('A','C','T','G')
        {
            let mut lmer_mutated = lmer.clone();
            lmer_mutated.insert(pos, c);
            //println!("b4 {} pos {} a8 {}",&lmer,pos, &lmer_mutated);

            assert_eq!(levenshtein(&lmer_mutated,&lmer), 1);
            res.insert(lmer_mutated.to_string());
        }
    }
    //println!("ball around {} is {:?}",lmer,res);
    res.iter().cloned().collect()
}

fn levenshtein_ball_2(lmer: &String) -> Vec<String>
{
    let ball_1 = levenshtein_ball_1(&lmer);
    let mut res = HashSet::new();
    for elt_from_ball_1 in ball_1
    {
        for elt_from_ball_2 in levenshtein_ball_1(&elt_from_ball_1)
        {
            if levenshtein(&elt_from_ball_2,&lmer) == 2
            { res.insert(elt_from_ball_2); }
        }
    }
    //println!("ball2 around {} is {:?}",lmer,res);
    res.iter().cloned().collect()
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

pub fn minhash_minimizer_decide(lmer: &str, params: &Params, lmer_counts: &mut HashMap<String,u32>) -> bool
{
    let size_miniverse = params.size_miniverse;
    let density = params.density;
    let l = params.l;
    let mut lmer_rev = String::new ();
    let mut lmerMin = String::new ();

    let h0 = city::hash32(&lmer) % size_miniverse;

    /* hash_mode == 1 {
       if start > 0
       {
       h1.remove(l, seq.as_bytes()[start-1]);
       h1.update(seq.as_bytes()[start + l-1]);
       }
       }*/
    
    if lmer_frequency_based
    {
        if revcomp_aware {
            lmer_rev = utils::revcomp(&lmer);
            lmerMin = std::cmp::min(lmer.to_string(), lmer_rev);
            println!("{}", lmerMin);
        }
        // minimizers must have high count
        let count = lmer_counts.entry(lmer.to_string()).or_insert(0);
        if (*count == 0) {
            return false;
        }
        let weight = (*count)/ (params.average_lmer_count as u32);
        return h0 < ((size_miniverse as f64*(density/(l as f64))) as u32 * weight);

        //if *count < params.average_lmer_count as u32 { return false;}
    }

    h0 < ((size_miniverse as f64*(density/(l as f64))) as u32 ) 
}

pub fn minhash(seq: &str, params: &Params, lmer_counts: &mut HashMap<String,u32>, minimizer_to_int: &HashMap<String,u32>) -> (Vec<String>, Vec<u32>, Vec<u32>)
{
    let l = params.l;
    
    //let hash_mode = 0; // 1: rolling (currently incompat with revcomp), 0: cityhash
    let mut res = Vec::new();
    let mut pos :Vec<u32> = Vec::new();
    if seq.len() < l { return (res, pos.clone(), pos.clone()); }

    //let mut h1 = RollingAdler32::from_buffer(&seq.as_bytes()[..l]);
    for start in 0..(&seq.len()-l+1)
    {
        let mut lmer = String::from(&seq[start..start+l]);
        if lmer.contains("N") { continue; }
        lmer = normalize_minimizer(&lmer);
        let mut should_insert : Option<String> = None;

        if params.levenshtein_minimizers == 0
        {
            if minhash_minimizer_decide(&lmer, params, lmer_counts)
            {
                should_insert = Some(lmer.to_string());
            }
        }

        // examine l-1, l, l+1 -mers, as a necessity in that scheme
        // small simplification: only one minimizer per start position
        if params.levenshtein_minimizers > 0
        {
            if minimizer_to_int.contains_key(&lmer) { should_insert = Some(lmer.to_string());}
            else
            {
                let lmer_minus_one = normalize_minimizer(&String::from(&seq[start..start+l-1]));
                if minimizer_to_int.contains_key(&lmer_minus_one) { should_insert = Some(lmer_minus_one.to_string());}
                else
                {
                    if start+l+1 < seq.len()
                    {
                        let lmer_plus_one = normalize_minimizer(&String::from(&seq[start..start+l+1]));
                        if minimizer_to_int.contains_key(&lmer_plus_one) { should_insert = Some(lmer_plus_one.to_string());}
                    }
                }
            }
        }

        if !should_insert.is_none()      
        {
            res.push(should_insert.unwrap());
            pos.push(start as u32);
            //println!("selected lmer");
        }
    }
        
    // convert minimizers to their integer representation
    let read_transformed : Vec<u32> = res.iter().map(|minim| minimizer_to_int[minim]).collect();

    (res, pos, read_transformed)
}

// TODO finish implementing it
#[allow(dead_code)]
fn wk_minimizers(seq: &str, l :usize) -> Vec<String>
{
    // actually it's (w,l)-minimizers, k is already used in the dbg
    let w = 15;
    let hash_mode = 0; // 0: cityhash
    let mut res = Vec::new();
    if seq.len() < l { return res; }
    let mut h0 :u32 = 0;
    let nl =  4f32.powf(l as f32) as u32;
    // TODO create a set of (minimizers,positions)

    for start in 0..(&seq.len()-w+1)
    {
        let mut min : u32 = std::u32::MAX;
        for pos in start..start+w
        {
            let mut lmer = String::from(&seq[start..start+l]);
            lmer = normalize_minimizer(&lmer);
            h0 = city::hash32(lmer) % nl;
            min = std::cmp::min(h0,min); // TODO record position as well
            // TODO add(minimizer,position) to set
        }

        //&res.push(lmer.to_string());
        res.push(h0.to_string()); // placeholder
        h0 = 1; // placeholder
    }

    // TODO convert the list of (minimizers,positions) to a vector of list minimizers
    res
}

// https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter

pub fn minimizers_preparation(lmer_counts: &mut HashMap<String,u32>, mut params: &mut Params, filename :&PathBuf, file_size: u64, levenshtein_minimizers: usize) -> (HashMap<String,u32>, HashMap<u32,String>) {

    let l = params.l;
    //let mut lmer_counts : HashMap<String,u32> = HashMap::new(); // for reference, 4^12 = 16M
    //if lmer_frequency_based {
    //    lmer_counting(&mut lmer_counts, &filename, file_size, &mut params);
    //}


    let mut list_minimizers : Vec<String> = Vec::new();
    
    // the following code replaces what i had before:
    // https://stackoverflow.com/questions/44139493/in-rust-what-is-the-proper-way-to-replicate-pythons-repeat-parameter-in-iter
    let multi_prod = (0..l).map(|i| vec!('A','C','T','G'))
            .multi_cartesian_product();

//    for lmer in kproduct("ACTG".to_string(), l as u32) {
    for lmer_vec in multi_prod {
        let lmer :String = lmer_vec.into_iter().collect();
        //println!("testing minimizer {}",lmer.to_string());
        if revcomp_aware {
            let lmer_rev = utils::revcomp(&lmer);
            if lmer > lmer_rev {continue;} // skip if not canonical
        }

        if ! minhash_minimizer_decide(&lmer, &params, lmer_counts) { 
            continue; 
        }
        //println!("found minimizer {}",lmer.to_string());

        list_minimizers.push(lmer);
    }
   
    let mut minimizer_to_int : HashMap<String,u32> = HashMap::new();
    let mut int_to_minimizer : HashMap<u32,String> = HashMap::new();
    let mut minim_idx : u32 = 0;

    if levenshtein_minimizers > 0
    {
        println!("adjusting for levenshtein minimizers");
        
        let mut set_nonoverlap_minimizers : HashSet<String> = HashSet::new();
        list_minimizers.sort(); // to make it deterministic
        for lmer in list_minimizers
        {
            let mut can_insert = true;
            let mut to_insert :Vec<String> = Vec::new();
        
            // otherwise consider all the lmers within his levenshtein-ball
            for lmer_neighbor_unmut in levenshtein_ball(&lmer, levenshtein_minimizers)
            {
                let lmer_neighbor = normalize_minimizer(&lmer_neighbor_unmut);

                if set_nonoverlap_minimizers.contains(&lmer_neighbor)
                {
                    can_insert = false;
                    break;
                }
                else
                {
                    to_insert.push(lmer_neighbor);
                }
            }

            if can_insert
            {
                for lmer_in_ball in to_insert
                {
                    // remember which minimizers we've already inserted
                    set_nonoverlap_minimizers.insert(lmer_in_ball.clone());

                    // assign numbers to minimizers, but each one in the ball has the same number
                    minimizer_to_int.insert(lmer_in_ball.to_string(),  minim_idx);
                }
                int_to_minimizer.insert(minim_idx,         lmer.to_string()); // bit of a hack
                minim_idx += 1;
            }
        }
    }
    else // regular good old minimizers
    {
        // assign numbers to minimizers, the regular way
        for lmer in list_minimizers
        {
            minimizer_to_int.insert(lmer.to_string(),  minim_idx);
            int_to_minimizer.insert(minim_idx,         lmer.to_string());
            minim_idx += 1;

            let count = lmer_counts.get(&lmer.to_string()).unwrap_or_else(|| &0);
            if *count != 0
            {
                //println!("found minimizer {} min-id {} count {}",lmer.to_string(),minim_idx, count);
            }
            else
            {
                //println!("found minimizer {} min-id {}",lmer.to_string(), minim_idx);
            }
        }
        assert_eq!(minimizer_to_int.len(), int_to_minimizer.len());
    }
    
    println!("selected {} minimizer ID's, {} sequences",int_to_minimizer.len(), minimizer_to_int.len());

    (minimizer_to_int, int_to_minimizer)
}

