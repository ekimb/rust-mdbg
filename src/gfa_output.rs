use petgraph::graph::DiGraph;
use petgraph::visit::EdgeRef;
use std::path::PathBuf;
use std::fs::File;
use std::error::Error;
use std::io::Write;
use super::Kmer;
use crate::kmer_vec::get;
use std::collections::HashMap;
use strsim::levenshtein;

use crate::utils::revcomp;

fn determine_orientation(id1: usize, id2: usize, nodes_vect: &Vec<&Kmer> ) -> (&'static str,&'static str)
{
    let n1 = nodes_vect[id1];
    let n2 = nodes_vect[id2];
    let rev_n1 = n1.reverse();
    let rev_n2 = n2.reverse();
    //println!("n1 {:?} n2 {:?} +/+? {}",n1.minimizers(),n2.minimizers(), n1.suffix() == n2.prefix());
    if n1.suffix() == n2.prefix() { return ("+","+"); }
    if n1.suffix() == rev_n2.prefix() { return ("+","-"); }
    if rev_n1.suffix() == n2.prefix() { return ("-","+"); }
    if rev_n1.suffix() == rev_n2.prefix() { return ("-","-"); }
    panic!("unknown orientation");
}   

fn matches(seq1 :&str, seq2 :&str, levenshtein_minimizers: usize) -> bool
{
    if levenshtein_minimizers == 0
    {
        return seq1.contains(seq2);
    }
    else
    {
        return levenshtein(seq1, seq2) <= levenshtein_minimizers+1 as usize; // +1 to take into account that some minimizers in the ball have length +- 1
    }
}

// should be fairly easy to find, as both sequences should share exactly k-1 minimizers exactly in
// common, and we even know which ones
fn find_overlap(seq1 :&str, seq2 :&str, ori1 :&str, ori2: &str, kmer1 :&Kmer, kmer2 :&Kmer, int_to_minimizer :&HashMap<u64,String>, minim_shift: &HashMap<Kmer,(usize,usize)>, levenshtein_minimizers: usize) -> u32
{
    // strategy: find the second minimizer at position 0 of seq2, by construction of the dbg
    /*
    let minim = kmer1.minimizers()[1];
    let minim_str = &int_to_minimizer[&minim];
    println!("looking for minim {} (rc: {})",  minim_str, revcomp(minim_str));
    let l = minim_str.len();

    let pos1_plus = seq1.find(minim_str);
    let pos1_minus = seq1.find(&revcomp(minim_str).to_string());

    let pos1 = match (pos1_plus,pos1_minus) {
        (None,None) => None,
        (Some(x),None) => Some((x,"+")),
        (None,Some(y)) => Some((y,"-")),
        (Some(x),Some(y)) => Some((x,"both")) 
    };
    assert!(!pos1.is_none());

    assert!(!pos2.is_none());

    let (pos1, ori1) = pos1.unwrap();
    let (pos2, ori2) = pos2.unwrap();

// problem with this strategy: the minimizer might appear at multiple positions in pos1

    println!("minim {}, seq1 {} seq2 {} min-pos1 {} min-ori1 {} min-pos2 {} min-ori2 {}",  minim_str, seq1, seq2, pos1,ori1,pos2,ori2);
    assert!(pos2 == 0); // sanity check
    (seq1.len() - pos1) as u32
    */
    
    // with minim_shift, it's relatively trivial now
    let shift_p :(usize,usize)= minim_shift[&kmer1.normalize().0];
    let shift_p :(usize,usize)= (shift_p.0 as usize, shift_p.1 as usize);
    //println!("minim shift {:?}",shift_p);
    /*for i in 0..kmer1.minimizers().len() {
        print!("{}\t", &int_to_minimizer[&kmer1.minimizers()[i]]);

    }
    print!("\n----------\t");
    for i in 0..kmer2.minimizers().len() {
        print!("{}\t", &int_to_minimizer[&kmer2.minimizers()[i]]);

    }*/
    let minim = kmer1.minimizers()[1];
    let minim_str = &int_to_minimizer[&minim];
    let l = minim_str.len() as usize;
    let minim_str_rev = revcomp(&minim_str);
    
    let mut shift :i32 = -1;
    //print!("seq {} minim {} (rc {}) ori1 {} shift_p {:?}\n",&seq1,minim_str,minim_str_rev,ori1,shift_p);
    if ori1 == "+" 
    {
        shift = shift_p.0 as i32; 
        if matches(&seq1[shift_p.0..shift_p.0+l], &minim_str, levenshtein_minimizers) || 
            matches(&seq1[shift_p.0..shift_p.0+l], &minim_str_rev, levenshtein_minimizers)
            { 
                shift = shift_p.0 as i32; 
                //print!("seq {} minim {} (rc {}) ori1 {} shift_p {:?}\n",&seq1,minim_str,minim_str_rev,ori1,shift_p);
            }
    }
    else
    {
        shift = shift_p.1 as i32;
        if matches(&seq1[shift_p.1..shift_p.1+l], &minim_str, levenshtein_minimizers) || 
            matches(&seq1[shift_p.1..shift_p.1+l], &minim_str_rev, levenshtein_minimizers)
            { 
                shift = shift_p.1 as i32;
                //print!("seq {} minim {} (rc {}) ori1 {} shift_p {:?}\n",&seq1,minim_str,minim_str_rev,ori1,shift_p);
            }
   }
    
    assert!(shift != -1);
    assert!((shift as usize) < seq1.len());
    shift as u32
}

pub fn output_gfa(gr: &DiGraph::<Kmer,Kmer>, dbg_nodes: &HashMap<Kmer,u32>, output_prefix :&PathBuf, kmer_seqs :&HashMap<Kmer,String>, int_to_minimizer :&HashMap<u64,String>, minim_shift: &HashMap<Kmer,(usize,usize)>, levenshtein_minimizers: usize, output_base_space: bool)  {
    // create a index->kmer index
    let nodes_vect : Vec<&Kmer> = dbg_nodes.keys().collect();
    let empty_str = String::new();
    
    let path = format!("{}{}",output_prefix.to_str().unwrap(),".gfa");
    let mut file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path, why.description()),
        Ok(file) => file,
    };
    //write!(file, "H\tVN:Z:2.0\n").expect("error writing GFA header");
     write!(file, "H\tVN:Z:1\n").expect("error writing GFA header");


    for node in gr.node_indices() {
        let idx = node.index();
        let seq = &kmer_seqs.get(nodes_vect[idx]).unwrap_or(&empty_str);
        let s_line = format!("S\t{}\t{}\tLN:i:{}\n",idx,seq,seq.len());
        write!(file, "{}", s_line).expect("error writing s_line");
       // let s_line = format!("S\t{}\t{}\t{}\n",idx,seq.len(),seq);
       // write!(file, "{}", s_line).expect("error writing s_line");
    }

    for e in gr.edge_references() {
        let id1 = e.source().index();
        let id2 = e.target().index();
        let (ori1, ori2) = determine_orientation(id1,id2,&nodes_vect);
        
        let mut kmer1 = nodes_vect[id1].clone();
        let mut kmer2 = nodes_vect[id2].clone();
        let mut seq1 = kmer_seqs.get(&kmer1).unwrap_or(&empty_str).clone();
        if ori1 == "-" {
            seq1 = revcomp(&seq1);
            kmer1 = kmer1.reverse();
        }
        let mut seq2 = kmer_seqs.get(&kmer2).unwrap_or(&empty_str).clone();
        if ori2 == "-" {
            seq2 = revcomp(&seq2);
            kmer2 = kmer2.reverse();
        }
        let mut shift = 0;
        if (output_base_space)
        {
            shift = find_overlap(&seq1, &seq2, ori1, ori2, &kmer1, &kmer2, int_to_minimizer, minim_shift, levenshtein_minimizers);
        }
        let mut overlap_length = seq1.len() as u32 -shift;

        //println!("seq1 len {} seq2 len {} overlap length {}", seq1.len(), seq2.len(), overlap_length);
        //if (overlap_length as usize) > seq2.len()
        //{
            //println!("kmer1 {:?}\nkmer2 {:?}\n", kmer1.minimizers(), kmer2.minimizers());
            //println!("seq1 {}\nseq2 {}\nori1 {} ori2 {}", seq1, seq2, ori1, ori2);
        //}
        //assert!((overlap_length as usize) < seq1.len() && (overlap_length as usize) < seq2.len());
        //
        // fun fact: overlap length may be slightly smaller than seq2.length actually
        // this is due to indels compensating the possible proximity of next minimizer
        // so for now we'll just do this dirty 'fix'
        // since anyway the importance of overlap length field in GFA is quite relative in this
        // pipeline
        overlap_length = std::cmp::min(overlap_length, (std::cmp::max(1,seq2.len())-1) as u32);
        let id1_beg = shift;
        let id1_end = shift+overlap_length;
        let id2_beg = 0;
        let id2_end = overlap_length;
        let mut id1_str = String::new();
        let mut id2_str = String::new();
        let mut ori_fin = "+";
        if ori1 != ori2 {ori_fin = "-";}
        if ori1 == "+" {id1_str = id1.to_string() + "+";} else { id1_str = id1.to_string() + "-";}

        //paf format 
//Query sequence name Query sequence length Query start (0-based) Query end (0-based) Relative strand: "+" or "-" Target sequence name Target sequence length Target start on original strand (0-based) Target end on original strand (0-based) Number of residue matches Alignment block length Mapping quality (0-255; 255 for missing)
        //let paf_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\n", id1_str, seq1.len(), id1_beg, id1_end, ori_fin, id2, seq2.len(), id2_beg, id2_end, overlap_length, overlap_length);
        
        
        
        //gfa2 : s1+ s2- b1 e1 b2 e2
        let l_line = format!("L\t{}\t{}\t{}\t{}\t{}M\n", id1, ori1, id2, ori2, overlap_length);
        write!(file, "{}", l_line).expect("error writing l_line");

        /*if overlap_length as usize == seq2.len()-1 {
            let e_line = format!("E\t*\t{}\t{}\t{}\t{}\t{}\t{}$\t*\n", id1_str, id2_str, id1_beg, id1_end, id2_beg, id2_end);
            write!(file, "{}", e_line).expect("error writing e_line");

        }
        else {
            let e_line = format!("E\t*\t{}\t{}\t{}\t{}$\t{}\t{}\t*\n", id1_str, id2_str, id1_beg, id1_end, id2_beg, id2_end);
            write!(file, "{}", e_line).expect("error writing e_line");
        }*/
        //write!(paf_file, "{}", paf_line).expect("error writing paf_line");

    }
}
