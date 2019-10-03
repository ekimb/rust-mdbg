//mod super::params;
use super::utils;
use super::Params;
use super::kproduct;
use super::revcomp_aware;
use super::minhash_minimizer_decide;
use super::levenshtein_minimizers; 
use std::collections::HashMap;
use std::collections::HashSet;
use pbr::ProgressBar;
use bio::io::fasta;
use std::path::PathBuf;
use strsim::levenshtein;

const preliminary_lmer_counting: bool = false;

fn lmer_counting(lmer_counts: &mut HashMap<String,u32>, filename :&PathBuf, file_size :u64, params: &mut Params) {
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

fn levenshtein_ball(lmer: &String) -> Vec<String>
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
        for c in vec!('A','C','T','G')
        {
            let mut lmer_mutated = lmer.clone();
            lmer_mutated.remove(pos);
            //println!("b4 {} pos {} a8 {}",&lmer,pos, &lmer_mutated);

            assert_eq!(levenshtein(&lmer_mutated,&lmer), 1);
            res.insert(lmer_mutated.to_string());
        }
    }
    //println!("ball around {} is {:?}",lmer,res);
    return res.iter().cloned().collect();
}



pub fn minimizers_preparation(mut params: &mut Params, filename :&PathBuf, file_size: u64) -> (HashMap<String,u32>, HashMap<u32,String>, HashMap<String,u32>) {

    let l = params.l;
    let mut lmer_counts : HashMap<String,u32> = HashMap::new(); // for reference, 4^12 = 16M
    if preliminary_lmer_counting {
        lmer_counting(&mut lmer_counts, &filename, file_size, &mut params);
    }


    let mut list_minimizers : Vec<String> = Vec::new();
    for lmer in kproduct("ACTG".to_string(), l as u32) {
        if revcomp_aware {
            let lmer_rev = utils::revcomp(&lmer);
            if lmer > lmer_rev {continue;} // skip if not canonical
        }

        if ! minhash_minimizer_decide(&lmer, &params, &lmer_counts) { 
            continue; 
        }

        list_minimizers.push(lmer);
        //let count = lmer_counts.get(&lmer.to_string()).unwrap_or_else(|| &0);
        //if *count != 0
        //{
            //println!("found minimizer {} count {}",lmer.to_string(),count);
        //}
        //println!("found minimizer {} min-id {}",lmer.to_string(), minim_idx);
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
        
            // otherwise consider all the lmers within his 1-levenshtein-ball
            for lmer_neighbor_unmut in levenshtein_ball(&lmer)
            {
                let mut lmer_neighbor = lmer_neighbor_unmut;
                if revcomp_aware {
                    let mut_rev = utils::revcomp(&lmer_neighbor);
                    lmer_neighbor = std::cmp::min(lmer_neighbor, mut_rev);
                }

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
    else
    {
        // assign numbers to minimizers, the regular way
        for lmer in list_minimizers
        {
            minimizer_to_int.insert(lmer.to_string(),  minim_idx);
            int_to_minimizer.insert(minim_idx,         lmer.to_string());
            minim_idx += 1;
        }
        assert_eq!(minimizer_to_int.len(), int_to_minimizer.len());
    }
    
    println!("selected {} minimizer ID's, {} sequences",int_to_minimizer.len(), minimizer_to_int.len());

    (minimizer_to_int, int_to_minimizer, lmer_counts)
}

