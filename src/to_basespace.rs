#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
#![allow(warnings)]
use pbr::ProgressBar;
use std::io::{self, stderr};
use std::error::Error;
use std::io::Write;
use std::io::{BufWriter, BufRead, BufReader};
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use itertools::Itertools;
use std::iter::FromIterator;
use std::fs::{File,remove_file};
use std::collections::HashSet;
extern crate array_tool;
//use adler32::RollingAdler32;
use std::fs;
use structopt::StructOpt;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use std::mem::{self, MaybeUninit};
use std::path::Path;
use lzzzz::lz4f::BufReadDecompressor;
use glob::glob;
mod utils;

#[derive(Debug, StructOpt)]
#[structopt(name = "to-basespace", about = "Converts the output of gfatools asm -u [any number of simplifications] to base-space.")]
struct Opt {
    /// Activate debug mode
    // short and long flags (-d, --debug) will be deduced from the field's name
    #[structopt(short, long)]
    debug: bool,
    /// Input file
    #[structopt(parse(from_os_str), short, long)]
    gfa: Option<PathBuf>,
    #[structopt(parse(from_os_str), short, long)]
    sequences: Option<PathBuf>
}

/// Try to get memory usage (resident set size) in bytes using the `getrusage()` function from libc.
// from https://github.com/digama0/mm0/blob/bebd670c5a77a1400913ebddec2c6248e76f90fe/mm0-rs/src/util.rs
fn get_memory_rusage() -> usize {
  let usage = unsafe {
    let mut usage = MaybeUninit::uninit();
    assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
    usage.assume_init()
  };
  usage.ru_maxrss as usize * 1024
}


// https://doc.rust-lang.org/rust-by-example/std_misc/file/read_lines.html
// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn read_lines_lz4<P>(filename: P) -> io::Result<io::Lines<BufReadDecompressor<'static, BufReader<File>>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(BufReadDecompressor::new(BufReader::new(file)).unwrap().lines())
}

fn main() {
    let start = Instant::now();
    let opt = Opt::from_args();      
    let mut gfa_file;
    let mut sequences_file;
    if !opt.gfa.is_none()       {gfa_file       = opt.gfa.unwrap();} 	   else {panic!("please specify an input GFA file (output of `gfatools asm [..] -u`).");} 
    if !opt.sequences.is_none() {sequences_file = opt.sequences.unwrap();} else {panic!("Please specify the prefix of [prefix].*.sequences files.");} 

    let debug = opt.debug;
    //let mut pb = ProgressBar::on(stderr(),file_size);
    // Step 1 : read the simplified GFA file
    //
    // we record which sequence we will need lin full and which we only need some extremity of
    let mut unitigs : HashMap<String, Vec<(u64, bool)>> = HashMap::new();
    let mut node2unitig : HashMap<u64, String> = HashMap::new();
    let mut unitig_name = String::new();
    let mut process_gfa_line = | line :&str, current_unitig : &mut (String, Vec<(u64, bool)>)| -> bool {
        //println!("line: {}", line);
        let is_S = line.starts_with('S');
        let is_L = line.starts_with('L');
        if is_S || is_L {
            //S       49194   *       LN:i:1  KC:i:157
            let mut v : Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
            unitig_name = v[1].clone();
            if current_unitig.1.len() > 0 {
                unitigs.insert(current_unitig.0.clone(), current_unitig.1.to_vec());
            }
            else {
                if is_L {return false;}
            }
            current_unitig.0 = unitig_name.clone(); // fill in the first unitig name
            current_unitig.1.clear();
        }
        else {
            if line.starts_with('A') {
                //A       utg0000001l     0       +       51909   0       1
                let v : Vec<&str> = line.split('\t').collect();
                let node_index = v[4].parse::<u64>().unwrap();
                current_unitig.1.push((node_index, v[3] == "+"));
                //println!("added node {} to node2unitig",node_index);
                node2unitig.insert(node_index, unitig_name.clone());
            }
        }
        return true;
    };
    
    if let Ok(lines) = read_lines(&gfa_file) {
	// Consumes the iterator, returns an (Optional) String
        let mut current_unitig : (String, Vec<(u64, bool)>) = ("".to_string(), Vec::new());
	    for line in lines {
	        if let Ok(line_contents) = line {
		        if !process_gfa_line(&line_contents, &mut current_unitig) {break;}
	        }
	    }
        // insert last unitig
        if current_unitig.1.len() > 0 {
            unitigs.insert(current_unitig.0.clone(), current_unitig.1.to_vec());
        }
    }
    println!("Done parsing GFA, got {} unitigs.", unitigs.len());

    // Step 1.5:
    // determine, for each node, whether to load its whole sequence or just the end

    enum LoadKind {
        Entire,
        EntireRc,
        Left,
        Right,
        LeftLast,
        RightLast
    };
    let mut load_node : HashMap<u64, LoadKind> = HashMap::new();
    for (unitig_name, unitig_vec) in &unitigs {
        for (i, (node_id, ori)) in unitig_vec.iter().enumerate() {
            if i == 0 {
                if *ori {load_node.insert(*node_id,LoadKind::Entire);} //+
                else {load_node.insert(*node_id,LoadKind::EntireRc);}
            }
            else {
                let is_last = i == unitig_vec.len() - 1;
                if *ori {load_node.insert(*node_id, if is_last {LoadKind::RightLast} else {LoadKind::Right});} // +
                else {load_node.insert(*node_id, if is_last {LoadKind::LeftLast} else {LoadKind::Left});}
            }
        }
    }
    //println!("Marked {} mDBG nodes to load sequences", load_node.len());
    // Step 2 : read the sequences file
    // we record in memory only the parts of the sequences we will need to fill in the complete GFA

    let mut sequences : HashMap<u64, String> = HashMap::new();
    let mut shifts : HashMap<String, ((u32, u32),(u32, u32))> = HashMap::new(); // also remember minim positions at the unitig-level for determining overlaps
    let mut process_sequence_line = |line: &str| {
        // 1740    [17721834324525059, 145164613009833572, 132671282722695729, 55654971230978398, 77488562567813534, 64734178350999831, 59815507667718894] TGGAAGCCTTTTTGATCTTCCTTACCACGGTGATGATGGTGTCGCTGATGATGGCCGCGGATCCCTCCCTGCTTGCTACGCCGCGTACGTACCTGATGAGCCATATGCCGTGGCTACCGTTTTTGCTGATCCTGCTGCCCGCCAACATCATGACCATGGTGATGTATGCCTTTCGTGCGGAACGCAAACACATTTCCGAAAGCGAAACCCACTTTCGGAACGCGATGGAATATTCCGCTATCGGTATGGCGTTAGTGGGCACCGAGGGACAATGGCTGCAAACCAACAAAGCGCTCTGCCAGT       *       *       (29, 5)
        if line.starts_with("#") {return;}
        let v : Vec<&str> = line.split('\t').collect();
        let node_id = v[0].parse::<u64>().unwrap();
        let unitig_name = node2unitig.get(&node_id);
        if !unitig_name.is_some() {return;}  // that node isn't used in a unitig.. weird.
        //println!("unitig {} seq {} minim {:?}",unitig_name.unwrap(),v[0],v[5][1..v[5].len()-1].split(',').collect::<Vec<&str>>());
        let minim_pos = v[5][1..v[5].len()-1].split(',').map(|s| s.trim().parse::<u32>().unwrap()).collect::<Vec<u32>>();
        let seq = v[2];
        let seq_len = seq.len() as u32;
        let left_seq = utils::revcomp(&v[2][0..minim_pos[0] as usize]);
        let right_seq = String::from(&v[2][seq.len()-minim_pos[1] as usize..]);
        let mut cur_shift = shifts.entry(unitig_name.unwrap().clone()).or_insert(((0, 0), (0, 0))).clone();
        match load_node.get(&node_id) {
            Some(load_kind) => {
                match load_kind {
                    LoadKind::Entire    =>    {sequences.insert(node_id, seq.to_string());         cur_shift.0 = (minim_pos[1], seq_len);}, // check those cur_shifts for correctness
                    LoadKind::EntireRc  =>    {sequences.insert(node_id, utils::revcomp(seq));     cur_shift.0 = (minim_pos[0], seq_len);},
                    LoadKind::Left      =>    {sequences.insert(node_id, left_seq.to_string());},
                    LoadKind::Right     =>    {sequences.insert(node_id, right_seq.to_string());},
                    LoadKind::LeftLast  =>    {sequences.insert(node_id, left_seq.to_string());    cur_shift.1 = (minim_pos[0], seq_len)},
                    LoadKind::RightLast =>    {sequences.insert(node_id, right_seq.to_string());   cur_shift.1 = (minim_pos[1], seq_len)},
                    _ => { panic!("Got unexpected load kind."); }
                }
            }
            _ => (),
        }
        shifts.insert(unitig_name.unwrap().clone(), cur_shift); // re-insert shifts in case they were modified
    };
    
    for path in glob(&format!("{}.*.sequences", &sequences_file.to_str().unwrap())).expect("Failed to read glob pattern for .sequences files.")  {
        let path = path.unwrap();
        let path = path.to_str().unwrap(); // rust really requires me to split the let statement in two..
        if let Ok(lines) = read_lines_lz4(path) {
            for line in lines {
                if let Ok(line_contents) = line {process_sequence_line(&line_contents);}
            }
        }
    }
    
    println!("Done parsing .sequences file, recorded {} sequences.", sequences.len());

    // Step 3 : read again the GFA file and at the same time write the .complete.gfa file with proper sequence

    let path = format!("{}{}", gfa_file.to_str().unwrap(),".complete.gfa");
    let mut complete_gfa_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}.", path, why.description()),
        Ok(file) => file,
    };
    write!(complete_gfa_file, "H\tVN:Z:1.0\n").expect("Error writing GFA header.");
    let reconstruct_seq = |unitig_name: &String| -> String {
        let mut res = String::new();
        //println!("reconstructing unitig {}",unitig_name);
        let unitig = unitigs.get(unitig_name).unwrap();
        for (node_id, ori) in unitig.iter() {
            let kminmer_seq = sequences.get(&node_id).unwrap();
            res.push_str(kminmer_seq);
        }
        res
    };
    let mut seq_lens : HashMap<String, usize> = HashMap::new();

    let mut process_gfa_line2 = |line: &str| {
        //println!("line: {}", line);
        let is_S = line.starts_with('S');
        let is_L = line.starts_with('L');
        let is_A = line.starts_with('A');
        if is_S {
            //S       49194   *       LN:i:1  KC:i:157
            //println!("{}",line);
            let mut v : Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
            let unitig_name = v[1].clone();
            let seq = reconstruct_seq(&unitig_name);
            v[2]= seq.clone();
            //v[3] = String::from(format!("LN:i:{}",seq.len())); // should already be there
            v[3] = String::from(format!("LN:i:{}", seq.len())); // but actually we want to fix it given that overlaps 
                                                               // were very approximately calculated. then gfatools complains
                                                               // and might even crash
            seq_lens.insert(unitig_name, seq.len());
            let s_line = v.join("\t");
            let s_line = v[..4].join("\t"); // deletes the RC/lc fields
            write!(complete_gfa_file, "{}\n", s_line).expect("Error writing S line.");
        }
        if is_L {
            //L       utg0000083l     -       utg0000084l     +       0M      L1:i:1  L2:i:1
            let mut v : Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
            // find the overlap length // actually it's maybe buggy. I ended up storing sequence
            // lengths and overlap lengths in the original GFA
            if false {
                let source_name = &v[1];
                let sink_name = &v[3];
                //println!("{}",source_name);
                let source_shift = shifts.get(source_name).unwrap();
                let source_shift = if v[2] == "+" {source_shift.0} else {source_shift.1}; // find_overlap() from complete_gfa.py
                let source_len = source_shift.1; //*seq_lens.get(source_name).unwrap() as u32;
                let source_shift = source_shift.0;
                let sink_shift = shifts.get(sink_name).unwrap();
                let sink_shift = if v[4] == "+" {sink_shift.0} else {sink_shift.1}; // find_overlap() from complete_gfa.py
                let sink_len   = sink_shift.1; //*seq_lens.get(sink_name).unwrap() as u32;
                let overlap_len =  std::cmp::min(source_len - source_shift, sink_len - 1);
                v[5] = String::from(format!("{}M", overlap_len));
            }
            let overlap_len = v[5][..v[5].len()-1].parse::<usize>().unwrap();
            let source_name = v[1].clone(); // am in a hurry
            let sink_name   = v[3].clone();
            if overlap_len > seq_lens[&source_name] || overlap_len > seq_lens[&sink_name] {
                let overlap_len = std::cmp::min(seq_lens[&source_name] - 1, seq_lens[&sink_name] - 1);
                // muted it because a bit too chatty, and not so useful
                //println!("fixing overlap for L-line {} (len:{}) {} (len:{}) from {} to {}M",source_name,seq_lens[&source_name],sink_name,seq_lens[&sink_name],v[5],overlap_len);
                v[5] = String::from(format!("{}M", overlap_len));
            }
            let l_line = v.join("\t");
            let l_line = v[..6].join("\t"); // deletes the L1:i and L2:i fields for now, I don't think they're essential
            write!(complete_gfa_file, "{}\n", l_line).expect("Error writing L line.");
        }
        if is_A {
            //A       utg0000001l     0       +       51909   0       1
            //let a_line = line; // copy as-is
            //write!(complete_gfa_file, "{}", a_line).expect("error writing A_line");
            // maybe let's not write it yet, until we compute the true offsets.
        }
        return true;
    };

    if let Ok(lines) = read_lines(&gfa_file) {
	// Consumes the iterator, returns an (Optional) String
	    for line in lines {
	        if let Ok(line_contents) = line {process_gfa_line2(&line_contents);}
	    }
    }
    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Maximum RSS: {:?} GB", (get_memory_rusage() as f32)/1024.0/1024.0/1024.0);
}

