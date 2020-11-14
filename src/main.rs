#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
#![allow(warnings)]
use pbr::ProgressBar;
use std::io::stderr;
use std::error::Error;
use std::io::Write;
use std::io::{BufWriter, BufRead, BufReader};
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use itertools::Itertools;
use closure::closure;
use std::iter::FromIterator;
use crate::kmer_vec::get;
use crate::read::Read;
use std::io::Read as OtherRead;
use std::fs::{File,remove_file};
use std::collections::HashSet;
extern crate array_tool;
use std::fs;
use crossbeam_utils::{thread};
use structopt::StructOpt;
use std::sync::{Arc, Mutex, MutexGuard};
use std::path::PathBuf;
use strsim::levenshtein;
use std::time::{Duration, Instant};
use std::mem::{self, MaybeUninit};
use editdistancewf as wf;
use seq_io::fasta;
use seq_io::core::BufReader as OtherBufReader;
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use lzzzz::lz4f::{WriteCompressor, ReadDecompressor, Preferences, PreferencesBuilder, CLEVEL_HIGH};
use xx_bloomfilter::Bloom;
use flate2::read::GzDecoder;
use std::sync::atomic::{AtomicUsize, Ordering};
use glob::glob;
use dashmap::{DashMap,mapref::one::Ref};
use thread_id;
use std::cell::UnsafeCell;
use std::io::Result;
mod utils;
mod gfa_output;
mod seq_output;
mod minimizers;
mod ec_reads;
mod kmer_vec;
mod poa;
mod read;
mod pairwise;
mod presimp;
use std::env;

const revcomp_aware: bool = true; // shouldn't be set to false except for strand-directed data or for debugging

type Kmer = kmer_vec::KmerVec;
type Overlap= kmer_vec::KmerVec;
type DbgIndex = u32;// heavily optimized assuming we won't get more than 2B kminmers of abundance <= 65535
type DbgAbundance = u16;
#[derive(Clone)] // TODO opt I only put this derive(Clone) to copy the DbgEntry in n2_entry in the code below, remove it when no more copies occur
struct DbgEntry { index: DbgIndex, abundance: DbgAbundance, seqlen: u32, shift: (u16,u16) } 
//type SeqFileType = WriteCompressor<File>;
struct SeqFileType (WriteCompressor<File>);
unsafe impl Sync for SeqFileType {} // same trick as below. we won't share files among threads but Rust can't know that.
impl SeqFileType {
    fn new(v: WriteCompressor<File>) -> SeqFileType   {
        SeqFileType(v)
    }
}
impl Write for SeqFileType {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        self.write(buf)
    }
    fn flush(&mut self) -> Result<()> {
        self.flush()
    }
}


struct RacyBloom(UnsafeCell<Bloom>); // intentionnally allowing data races as a tradeoff for bloom speed
unsafe impl Sync for RacyBloom {}

// follows https://sodocumentation.net/rust/topic/6018/unsafe-guidelines
// don't try this at home
impl RacyBloom {
    fn new(v: Bloom) -> RacyBloom {
        RacyBloom(UnsafeCell::new(v))
    }
    fn get(&self) -> &mut Bloom {
        // UnsafeCell::get() returns a raw pointer to the value it contains
        // Dereferencing a raw pointer is also "unsafe"
        unsafe { &mut *self.0.get() }
    }
}

type ThreadIdType = usize;

pub struct Params
{
    l: usize,
    k: usize,
    n: usize,
    t: usize,
    w: usize,
    density :f64,
    size_miniverse: u32,
    average_lmer_count : f64,
    lmer_counts_min : u32,
    lmer_counts_max : u32,
    min_kmer_abundance : DbgAbundance,
    levenshtein_minimizers : usize,
    correction_threshold: i32,
    distance: usize ,
    reference: bool,
    uhs: bool,
    error_correct: bool,
    has_lmer_counts: bool,
    debug: bool,
}


fn debug_output_read_minimizers(seq_str: &String, read_minimizers : &Vec<String>, read_minimizers_pos :&Vec<u32>)
{
    println!("\nseq: {}",seq_str);
    print!("min: ");
    let mut current_minimizer :String = "".to_string();
    for i in 0..seq_str.len()
    {
        if read_minimizers_pos.contains(&(i as u32))
        {
            let index = read_minimizers_pos.iter().position(|&r| r == i as u32).unwrap();
            current_minimizer = read_minimizers[index].clone();
            let c = current_minimizer.remove(0);
            if c == seq_str.chars().nth(i).unwrap()
            {
                print!("X");
            }
                else
            {
                print!("x");
            }
            continue;
        }
        if current_minimizer.len() > 0
        {
            let c = current_minimizer.remove(0);
            print!("{}",c);
        }
        else
        {
            print!(".");
        }
    }
    println!("");

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


// thread helpers
fn thread_update_hashmap<U,V>(hashmap_all: &Arc<Mutex<HashMap<usize,HashMap<U,V>>>>, hashmap: HashMap<U,V>, thread_num: usize)
{
    let mut hashmap_all = hashmap_all.lock().unwrap();
    let mut entry = hashmap_all.entry(thread_num).or_insert(HashMap::new());
    *entry = hashmap; // I believe hashmap is moved in this function as per https://stackoverflow.com/a/29490907 
}

pub fn thread_update_vec<U>(vec_all: &Arc<Mutex<HashMap<usize,Vec<U>>>>, vec: Vec<U>, thread_num: usize)
{
    let mut vec_all = vec_all.lock().unwrap();
    let mut entry = vec_all.entry(thread_num).or_insert(Vec::new());
    *entry = vec;
}

fn get_reader(path: &PathBuf) -> Box<OtherRead + Send> {
    let filetype = "gz";
    let file = match File::open(path) {
            Ok(file) => file,
            Err(error) => panic!("There was a problem opening the file: {:?}", error),
        };

    let reader: Box<OtherRead + Send> = match filetype { 
        "gz" => Box::new(GzDecoder::new(file)) as Box<dyn OtherRead + Send>, 
        _ => Box::new(file) as Box<dyn OtherRead + Send>, 
        
    }; 
    reader
}



#[derive(Debug, StructOpt)]
#[structopt(name = "rust-mhdbg", about = "Original implementation of MinHash de Bruijn graphs")]
struct Opt {
    /// Activate debug mode
    // short and long flags (-d, --debug) will be deduced from the field's name
    #[structopt(short, long)]
    debug: bool,

    /// Input file
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,
    #[structopt(long)]
    uhs: Option<String>,

    /// Output graph/sequences prefix 
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,

    #[structopt(short, long)]
    k: Option<usize>,
    #[structopt(short, long)]
    l: Option<usize>,
    #[structopt(short, long)]
    n: Option<usize>,
    #[structopt(short, long)]
    t: Option<usize>,
    #[structopt(long)]
    density: Option<f64>,
    #[structopt(long)]
    minabund: Option<usize>,
    #[structopt(short, long)]
    w: Option<usize>,
    #[structopt(long)]
    distance: Option<usize>,
    #[structopt(long)]
    correction_threshold: Option<i32>,
    #[structopt(long)]
    levenshtein_minimizers: Option<usize>,
    #[structopt(long)]
    test1: bool,
    #[structopt(long)]
    test2: bool,
    #[structopt(long)]
    no_error_correct: bool,
    #[structopt(long)]
    restart_from_postcor: bool,
    #[structopt(long)]
    reference: bool,
    #[structopt(parse(from_os_str), long)]
    lmer_counts: Option<PathBuf>,
    #[structopt(long)]
    lmer_counts_min: Option<u32>,
    #[structopt(long)]
    lmer_counts_max: Option<u32>,
    #[structopt(long)]
    presimp: Option<f32>,
    #[structopt(long)]
    threads: Option<usize>,

}


fn main() {
    let start = Instant::now();
    let opt = Opt::from_args();      
    let mut uhs : bool = false;
    let mut filename = PathBuf::new();
    let mut lmer_counts_filename = PathBuf::new();
    let mut uhs_filename = String::new();
    let mut output_prefix;
    let mut k: usize = 10;
    let mut l: usize = 12;
    let mut n: usize = 2;
    let mut t: usize = 0;
    let mut w: usize = 0;
    let mut density :f64 = 0.10;
    let mut min_kmer_abundance: DbgAbundance = 2;
    let mut levenshtein_minimizers: usize = 0;
    let mut distance: usize = 0;
    let mut error_correct: bool = true;
    let mut restart_from_postcor: bool = false;
    let mut correction_threshold : i32 = 0;
    let mut reference : bool = false;
    let mut windowed : bool = false;
    let mut has_lmer_counts : bool = false;
    let mut lmer_counts_min: u32 = 2;
    let mut lmer_counts_max: u32 = 100000;
    let mut presimp :f32 = 0.0;
    let mut threads : usize = 8;
    if opt.no_error_correct {
        error_correct = false;
    }
    if opt.reference {
        reference = true;
        error_correct = false; // rayan->baris: remove this comment once you see it. added this line to make sure we don't do POA on reference
    }
    if opt.test1 {
        filename = PathBuf::from("../read50x_ref10K_e001.fa"); 
        if opt.k.is_none() { k = 5; }
        if opt.l.is_none() { l = 8; }
        if opt.density.is_none() { density = 1.0 };
    }
    else if opt.test2 {
        filename = PathBuf::from("../SRR9969842_vs_chr4.fasta");
        if opt.k.is_none() { k = 50; }
        if opt.l.is_none() { l = 12; }
        if opt.density.is_none() { density = 0.1 };
    }
    if !opt.k.is_none() { k = opt.k.unwrap() } else { println!("Warning: using default k value ({})",k); } 
    if !opt.l.is_none() { l = opt.l.unwrap() } else { println!("Warning: using default l value ({})",l); }
    if !opt.n.is_none() { n = opt.n.unwrap() } else { println!("Warning: using default n value ({})",n); }
    if !opt.t.is_none() { t = opt.t.unwrap() } else { println!("Warning: using default t value ({})",t); }

    if !opt.density.is_none() { density = opt.density.unwrap() } else { println!("Warning: using default minhash density ({}%)",density*100.0); }
    if !opt.minabund.is_none() { min_kmer_abundance = opt.minabund.unwrap() as DbgAbundance } else { println!("Warning: using default min kmer abundance value ({})",min_kmer_abundance); }
    if !opt.w.is_none() { windowed = true; w = opt.w.unwrap(); } else { println!("Warning: Using default density-based"); }
    if !opt.presimp.is_none() { presimp = opt.presimp.unwrap(); } else { println!("Warning: Using default no-presimp"); }
    if !opt.threads.is_none() { threads = opt.threads.unwrap(); } else { println!("Warning: Using default num threads (8)"); }
    if !opt.correction_threshold.is_none() { correction_threshold = opt.correction_threshold.unwrap() } else { println!("Warning: using default correction threshold value ({})",correction_threshold); }

    if !opt.levenshtein_minimizers.is_none() { levenshtein_minimizers = opt.levenshtein_minimizers.unwrap() }
    if !opt.distance.is_none() { distance = opt.distance.unwrap() }
    if distance > 2 {distance = 2;}
    let distance_type = match distance { 0 => "jaccard", 1 => "containment", 2 => "mash",_ => "mash" };
    let minimizer_type = match levenshtein_minimizers { 0 => "reg", 1 => "lev1", 2 => "lev2",_ => "levX" };
    if opt.levenshtein_minimizers.is_none() { println!("Warning: using default minimizer type ({})",minimizer_type); }
    if opt.distance.is_none() { println!("Warning: using default distance metric ({})",distance_type); }
    if opt.restart_from_postcor { restart_from_postcor = true;}


    output_prefix = PathBuf::from(format!("{}graph-k{}-p{}-l{}",minimizer_type,k,density,l));

    if !opt.reads.is_none() { filename = opt.reads.unwrap().clone(); } 
    if !opt.lmer_counts.is_none() { 
        has_lmer_counts = true;
        lmer_counts_filename = opt.lmer_counts.unwrap().clone(); 
        if !opt.lmer_counts_min.is_none() { lmer_counts_min = opt.lmer_counts_min.unwrap(); } else { println!("Warning: Using default l-mer min count ({})", lmer_counts_min); }
        if !opt.lmer_counts_max.is_none() { lmer_counts_max = opt.lmer_counts_max.unwrap(); } else { println!("Warning: Using default l-mer max count ({})", lmer_counts_max); }
    } 

    if !opt.uhs.is_none() { 
        uhs = true;
        uhs_filename = opt.uhs.unwrap(); 
    } 
    if !opt.prefix.is_none() { output_prefix = opt.prefix.unwrap(); } else { println!("Warning: using default prefix ({})",output_prefix.to_str().unwrap()); }
    
    if filename.as_os_str().is_empty() { panic!("please specify an input file"); }

    let debug = opt.debug;

    let size_miniverse = match revcomp_aware
    {
        false => 4f32.powf(l as f32) as u32,
        true => 4f32.powf(l as f32) as u32 / 2
    };

    let mut params = Params { 
        l,
        k,
        n,
        t,
        w,
        density,
        size_miniverse,
        average_lmer_count: 0.0,
        lmer_counts_min,
        lmer_counts_max,
        min_kmer_abundance,
        levenshtein_minimizers,
        distance,
        correction_threshold,
        reference,
        uhs,
        error_correct,
        has_lmer_counts,
        debug,
    };
    // init some useful objects
    let mut nb_minimizers_per_read : f64 = 0.0;
    let mut nb_reads : u64 = 0;
    // get file size for progress bar
    let metadata = fs::metadata(&filename).expect("error opening input file");
    let file_size = metadata.len();
    let mut pb = ProgressBar::on(stderr(),file_size);
    let mut lmer_counts : HashMap<String, u32> = HashMap::new();

    if has_lmer_counts {
        let lmer_counts_file = match File::open(lmer_counts_filename) {
            Err(why) => panic!("couldn't load l-mer counts file: {}", why.description()),
            Ok(lmer_counts_file) => lmer_counts_file,
        }; 
        let mut br = BufReader::new(lmer_counts_file);
        loop
        {
            let mut line = String::new();
            let new_line = |line: &mut String, br :&mut BufReader<File>| { line.clear(); br.read_line(line).ok(); };
            if let Err(e) = br.read_line(&mut line) { break; }
            if line.len() == 0                      { break; }
            let trimmed  = line.trim().to_string();   
            let vec : Vec<String> = trimmed.split(" ").map(String::from).collect();
            let lmer = vec[0].to_string();
            let lmer_rev = utils::revcomp(&lmer);
            let lmer = if lmer > lmer_rev {lmer} else {lmer_rev}; //don't trust the kmer counter to normalize like we do
            let count = vec[1].parse::<u32>().unwrap();
            lmer_counts.insert(lmer, count);               
            new_line(&mut line, &mut br);
        }
    }
    
    // did some l-mer counting before, but not using this anymore.
    // the above code uses already-counted lmers
    //let mut lmer_counts : HashMap<String, u32> = HashMap::new();
    //minimizers::lmer_counting(&mut lmer_counts, &filename, file_size, &mut params);
 

    let mut minimizer_to_int: HashMap<String,u64> = HashMap::new();
    let mut int_to_minimizer: HashMap<u64,String> = HashMap::new();


    // only need to initialize the minimizer_to_int / int_to_minimizer array if we do POA or use robust minimizers
    // they can be costly for k=14
    if has_lmer_counts || error_correct {
        let res = minimizers::minimizers_preparation(&mut params, &filename, file_size, levenshtein_minimizers, &lmer_counts);
        minimizer_to_int = res.0;
        int_to_minimizer = res.1;
    }
    let mut fasta_reads : bool = false;
    
    let queue_len = 200; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html

    let mut uhs_kmers = HashMap::<String, u32>::new();
    if params.uhs {
        uhs_kmers = minimizers::uhs_preparation(&mut params, &uhs_filename)
    }

    // dbg_nodes is a hash table containing (kmers -> (index,count))
    // it will keep only those with count > 1
    let mut dbg_nodes     : Arc<DashMap<Kmer,DbgEntry>> = Arc::new(DashMap::new()); // it's a Counter
    let mut bloom : RacyBloom = RacyBloom::new(Bloom::new_with_rate(100_000_000, 1e-6)); // a bf to avoid putting stuff into dbg_nodes too early
    static node_index: AtomicUsize = AtomicUsize::new(0); // associates a unique integer to each dbg node
    let mut kmer_seqs     : HashMap<Kmer,String> = HashMap::new(); // associate a dBG node (k-min-mer) to an arbitrary sequence from the reads
    let mut kmer_seqs_lens: HashMap<Kmer,Vec<u32>> = HashMap::new(); // associate a dBG node to the lengths of all its sequences from the reads
    let mut kmer_origin   : HashMap<Kmer,String> = HashMap::new(); // remember where in the read/refgenome the kmer comes from, for debugging only

    // correction stuff
    let mut buckets : HashMap<Vec<u64>, Vec<String>> = HashMap::new();
    let mut corrected_map = HashMap::<String, (String, Vec<String>, Vec<usize>, Vec<u64>)>::new(); // reduce runtime of POA by simultaneous correction of template and its aligned reads. keeps track of reads that have been already corrected
    let postcor_path = PathBuf::from(format!("{}.postcor",output_prefix.to_str().unwrap()));
    let poa_path     = PathBuf::from(format!("{}.poa",    output_prefix.to_str().unwrap()));
    let mut reads_by_id = HashMap::<String, Read>::new();

    // delete all previous sequence files
    for path in glob(&format!("{}*.sequences", output_prefix.to_str().unwrap()).to_string()).expect("Failed to read glob pattern")  {
        let path = path.unwrap();
        let path = path.to_str().unwrap(); // rust really requires me to split the let statement in two..
        println!("removing old sequences file: {}",&path);
        fs::remove_file(path);
    }

    let mut seq_write = |file : &mut SeqFileType,s|
    {
        write!(file,"{}",s);
    };

    let mut sequences_files : Arc<DashMap<ThreadIdType,SeqFileType>> = Arc::new(DashMap::new());
    let create_sequences_file = |thread_id:ThreadIdType| -> SeqFileType {
        let seq_path = PathBuf::from(format!("{}.{}.sequences", output_prefix.to_str().unwrap(), thread_id));
        let mut file = match File::create(&seq_path) {
            Err(why) => panic!("couldn't create file: {}", why.description()),
            Ok(file) => file,
        };
        //let mut sequences_file = BufWriter::new(file);
        let mut sequences_file = SeqFileType::new(WriteCompressor::new(file, Preferences::default()).unwrap()); // regular lz4f
        //let mut sequences_file = WriteCompressor::new(&mut file, PreferencesBuilder::new().compression_level(CLEVEL_HIGH).build()).unwrap();  // too slow
        seq_write(&mut sequences_file,format!("# k = {}\n",k));
        seq_write(&mut sequences_file,format!("# l = {}\n",l));
        seq_write(&mut sequences_file,"# structure of remaining of the file:\n".to_string());
        seq_write(&mut sequences_file,"# [node name]\t[list of minimizers]\t[sequence of node]\t[abundance]\t[origin]\t[shift]\n".to_string());
        sequences_file
    };

    let mut add_kminmers = |vec: &Vec<(Kmer,String,bool,String,(usize,usize))>| 
    {
        let thread_id :usize =  thread_id::get();
        // determine to which sequence files to write
        if sequences_files.get(&thread_id).is_none()  {
            sequences_files.insert(thread_id, create_sequences_file(thread_id));
        }
        let sequences_file = sequences_files.get(&thread_id).unwrap().value();

        for (node, seq, seq_reversed, origin, shift) in vec.iter() {
            let mut abundance: u16 = 0;
            // 1) check if the kminmer is in the BF
            unsafe {
                if ! bloom.get().check_and_add(&node)   { 
                    // if not, necessarily it wasn't seen before so abundance is 1, and insert in BF
                    // and then, if we're dealing with a reference genome, insert in dbg_nodes
                    abundance = 1;
                } else {
                    // it was in the BF. get its true abundance
                    abundance = 2; // assume since it's in the BF that the abundance is >= 2 regardless of whether it's in dbg_nodes
                    // but maybe it's already in dbg_nodes, let's check
                    // TODO optimize that code
                    if dbg_nodes.contains_key(&node)  {
                        abundance = dbg_nodes.get(&node).unwrap().abundance;
                    }
                }
            }
            //println!("abundance: {}",abundance);
            if (params.reference && abundance == 1) || ((!params.reference) && abundance == min_kmer_abundance) {
                // only save each kminmer once, once we're sure it'll pass the filter.
                // (note that this doesnt enable to control which seq we save based on
                // median seq length)
                let seq = if *seq_reversed { utils::revcomp(seq) } else { seq.to_string() };
                let cur_node_index :DbgIndex = node_index.fetch_add(1,Ordering::Relaxed) as DbgIndex;
                let seq_line = format!("{}\t{}\t{}\t{}\t{}\t{:?}",cur_node_index,node.print_as_string(), seq, "*", origin, shift);

                seq_write(&mut sequences_files.get_mut(&thread_id).unwrap(), format!("{}\n", seq_line));

                // now record what we just saved
                let lowprec_shift = (shift.0 as u16, shift.1 as u16);
                dbg_nodes.insert(node.clone(),DbgEntry{index: cur_node_index, abundance: 2, seqlen: seq.len() as u32, shift: lowprec_shift}); 
            }
            else  {
                // at this point the element _has_ to already have been inserted in dbg_nodes
                // just increase its abundance 
                dbg_nodes.get_mut(&node).unwrap().abundance += 1; // if that code crashes it has to be a race condition i didnt check carefully
            }
        }
    };


    if ! restart_from_postcor
    {
        let mut ec_file  = BufWriter::new(File::create("/dev/null").unwrap());
        if error_correct || reference { // write this file even for an uncorrected reference, needed by evaluate_ec.py
            ec_file = ec_reads::new_file(&output_prefix); // reads before correction
        }
        
        // worker thread
        let process_read_aux = |seq_str: &str, seq_id: &str | -> (Vec<(Kmer,String,bool,String,(usize,usize))>,Read) {
            let mut output : Vec<(Kmer,String,bool,String,(usize,usize))> = Vec::new();
            let mut read_obj = Read::extract(&seq_id, &seq_str, &params, &minimizer_to_int, &int_to_minimizer);
            println!("Received read in worker thread, transformed len {}", read_obj.transformed.len());
            if read_obj.transformed.len() > k {
                output = read_obj.read_to_kmers(&params);
            }

            add_kminmers(&output); 

            (output, read_obj)
        };
        let process_read_fasta = |record: seq_io::fasta::RefRecord, found : &mut (Vec<(Kmer,String,bool,String,(usize,usize))>,Read) | {
            let seq_str = String::from_utf8_lossy(record.seq()).to_string(); // might induce a copy? can probably be optimized (see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html)
            let seq_id = record.id().unwrap().to_string();
            *found = process_read_aux(&seq_str,&seq_id);
        };
        let process_read_fastq = |record: seq_io::fastq::RefRecord, found : &mut (Vec<(Kmer,String,bool,String,(usize,usize))>,Read) | {
            let seq_str = String::from_utf8_lossy(record.seq()).to_string(); // might induce a copy? can probably be optimized (see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html)
            let seq_id = record.id().unwrap().to_string();
            *found = process_read_aux(&seq_str,&seq_id);
        };


        // parallel fasta parsing, with a main thread that writes to disk and populates hash tables
        let mut main_thread =  |found: &(Vec<(Kmer,String,bool,String,(usize,usize))>,Read)| { // runs in main thread
            nb_reads += 1;
            let (vec, read_obj) = found;
            println!("Received read in main thread, nb kmers: {}", vec.len());

            if error_correct || reference
            {
                reads_by_id.insert(read_obj.id.to_string(), read_obj.clone());

                ec_reads::record(&mut ec_file, &read_obj.id.to_string(), &read_obj.seq, &read_obj.transformed.to_vec(), &read_obj.minimizers, &read_obj.minimizers_pos);

                for i in 0..read_obj.transformed.len()-n+1 {
                    let n_mer = utils::normalize_vec(&read_obj.transformed[i..i+n].to_vec());
                    let mut entry = buckets.entry(n_mer).or_insert(Vec::<String>::new());
                    entry.push(read_obj.id.to_string());
                }
            }

            // Some(value) will stop the reader, and the value will be returned.
            // In the case of never stopping, we need to give the compiler a hint about the
            // type parameter, thus the special 'turbofish' notation is needed.
            None::<()>
        };

        let filename_str = filename.to_str().unwrap();
        if filename_str.ends_with(".fasta.gz") || filename_str.ends_with(".fa.gz")  || filename_str.ends_with(".fa") || filename_str.ends_with(".fasta") {
            fasta_reads = true;
            println!("Input file: {}", filename_str);
            println!("Format: FASTA");
        }
        let buf = BufReader::new(get_reader(&filename));
        println!("Parsing input sequences...");

        if fasta_reads {
            let reader = seq_io::fasta::Reader::new(buf);
            read_process_fasta_records(reader, threads as u32, queue_len, process_read_fasta, |record,found| { main_thread(found) });
        }
        else {
            let reader = seq_io::fastq::Reader::new(buf);
            read_process_fastq_records(reader, threads as u32, queue_len, process_read_fastq, |record,found| { main_thread(found) });
        }

        pb.finish_print("Done converting reads to k-min-mers.");
        println!("Number of reads: {}", nb_reads);


        // this part will correct the reads, and dump them to disk
        if error_correct {
            ec_reads::flush(&mut ec_file);
            let chunks = ec_reads::load(&output_prefix);
            let mut chunk_length = 1;
            let mut ec_entries      = Arc::new(Mutex::new(HashMap::<usize, Vec<(String, String, Vec<u64>, Vec<String>, Vec<usize>)>>::new()));
            let mut poa_entries     = Arc::new(Mutex::new(HashMap::<usize, HashMap<String, Vec<String>>>::new()));

            let mut ec_file_poa     = ec_reads::new_file(&poa_path);      // POA debug info (which reads were recruited per template, I think. Baris can correct/confirm)
            let mut ec_file_postcor = ec_reads::new_file(&postcor_path);  // reads after correction

            if chunks.len() > threads {chunk_length = chunks.len()/threads+1;}
            thread::scope(|s| {
                let mut guards = Vec::with_capacity(threads);
                for (thread_num, chunk) in chunks.chunks(chunk_length).enumerate() {
                    let mut ec_entries = ec_entries.clone();
                    let mut poa_entries = poa_entries.clone();

                    let guard = s.spawn(closure!(move chunk, ref params, ref int_to_minimizer, ref output_prefix, ref mut pb, ref buckets, ref reads_by_id, |_| {
                        let mut ec_entry = Vec::<(String, String, Vec<u64>, Vec<String>, Vec<usize>)>::new();
                        let mut corrected_map = HashMap::new();
                        let mut poa_map = HashMap::new();
                        for ec_record in chunk.iter() {
                            let mut read_obj = Read {id: ec_record.seq_id.to_string(), minimizers: ec_record.read_minimizers.to_vec(), minimizers_pos: ec_record.read_minimizers_pos.to_vec(), transformed: ec_record.read_transformed.to_vec(), seq: ec_record.seq_str.to_string(), corrected: false};
                            if !corrected_map.contains_key(&read_obj.id) { 
                                read_obj.poa_correct(&int_to_minimizer, &mut poa_map, &buckets, &params, &mut corrected_map, &reads_by_id);
                            }
                            else {
                                continue;
                                read_obj.seq = corrected_map[&read_obj.id].0.to_string();
                                read_obj.minimizers = corrected_map[&read_obj.id].1.to_vec();
                                read_obj.minimizers_pos = corrected_map[&read_obj.id].2.to_vec();
                                read_obj.transformed = corrected_map[&read_obj.id].3.to_vec();
                            }
                            ec_entry.push((read_obj.id.to_string(), read_obj.seq.to_string(), read_obj.transformed.to_vec(), read_obj.minimizers.to_vec(), read_obj.minimizers_pos.to_vec()));
                        }

                        thread_update_vec(    &ec_entries,  ec_entry, thread_num);
                        thread_update_hashmap(&poa_entries, poa_map, thread_num);

                    }));
                    guards.push(guard);
                }
            }).unwrap();

            pb.finish_print("Done with correction.");
            let mut ec_entries = ec_entries.lock().unwrap();
            let mut poa_entries = poa_entries.lock().unwrap(); 
            for thread_num in 0..threads {
                let mut ec = ec_entries.entry(thread_num).or_insert(Vec::new());
                for tuple in ec.iter() {
                    ec_reads::record(&mut ec_file_postcor, &tuple.0, &tuple.1, &tuple.2, &tuple.3, &tuple.4);
                }
                let mut poa = poa_entries.entry(thread_num).or_insert(HashMap::new());
                for (temp, vec) in poa.into_iter() {ec_reads::record_poa(&mut ec_file_poa, temp, vec.to_vec());}
            }
            ec_reads::flush(&mut ec_file_postcor);
            ec_reads::flush(&mut ec_file_poa);
        }
    }

    // now load the error corrected reads and construct the dbg from it
    // it's a little slow because read_to_kmers is called single threaded, there's room for
    // introducing multithreading here
    if error_correct || restart_from_postcor
    {
        dbg_nodes.clear(); // will be populated by add_kminmers()
        let chunks = ec_reads::load(&postcor_path);

        let mut sequences_file = create_sequences_file(0);

        node_index.store(0,Ordering::Relaxed);
        for ec_record in chunks.iter() {
            let mut read_obj = Read {id: ec_record.seq_id.to_string(), minimizers: ec_record.read_minimizers.to_vec(), minimizers_pos: ec_record.read_minimizers_pos.to_vec(), transformed: ec_record.read_transformed.to_vec(), seq: ec_record.seq_str.to_string(), corrected: false};
            if read_obj.transformed.len() > k { 
                let output = read_obj.read_to_kmers(&params);
                add_kminmers(&output);
            }
        }
        sequences_file.flush().unwrap();
    }

    for mut sequences_file in sequences_files.iter_mut() {
        sequences_file.flush().unwrap();
    }

    // now DBG creation can start
    println!("nodes before abund-filter: {}", dbg_nodes.len());
    dbg_nodes.retain(|x,c| c.abundance >= (min_kmer_abundance as DbgAbundance));
    println!("nodes after: {}", dbg_nodes.len());
    let mut dbg_edges : Vec<(&Kmer,&Kmer)> = Vec::new();

    let path = format!("{}{}",output_prefix.to_str().unwrap(),".gfa");
    let mut gfa_file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path, why.description()),
        Ok(file) => file,
    };
    write!(gfa_file, "H\tVN:Z:1\n").expect("error writing GFA header");

    // index k-1-mers
    // TODO opt: used to be a Vec<&Kmer> but unfortunately dashmap
    // doesnt seem to let me keep references of its keys.
    let mut km_index : HashMap<Overlap,Vec<Kmer>> = HashMap::new(); 
    for item in dbg_nodes.iter()
    {
        let node = item.key();
        let entry = item.value();
        // take this iteration opportunity to write S lines
        let length = entry.seqlen;
        let s_line = format!("S\t{}\t{}\tLN:i:{}\tKC:i:{}\n",entry.index,"*",length,entry.abundance);
        write!(gfa_file, "{}", s_line).expect("error writing s_line");

        let first :Overlap = node.prefix().normalize().0;
        let second  :Overlap = node.suffix().normalize().0;
        let mut insert_km = |key,val| 
        {
            match km_index.entry(key) {
                Entry::Vacant(e) => { e.insert(vec![val]); },
                Entry::Occupied(mut e) => { e.get_mut().push(val); }
            }
        };
        insert_km(first,node.clone());
        insert_km(second,node.clone());
    }

    let mut nb_edges = 0;
    let mut presimp_removed = 0;
    let mut removed_edges : HashSet<(DbgIndex,DbgIndex)> = HashSet::new();
    let mut vec_edges = Vec::new(); // for presimp, keep a list of all edges to insert

    // create a vector of dbg edges (if we want to create a petgraph)
    // otherwise just output edges to gfa without keeping them in mem
    for item in dbg_nodes.iter() 
    {
        let n1 = item.key();
        let n1_entry = item.value();
        let rev_n1 = n1.reverse();
        let n1_abundance = n1_entry.abundance;
        let n1_index     = n1_entry.index;
        let n1_seqlen    = n1_entry.seqlen;

        // bit of a rust noob way to code this, because i'm not too familiar with types yet..
        let key1=n1.suffix().normalize().0;
        let key2=n1.prefix().normalize().0;
        for key in [key1,key2].iter()
        {
            if km_index.contains_key(&key)
            {
                let list_of_n2s : &Vec<Kmer> = km_index.get(&key).unwrap();
                let mut potential_edges : Vec<(DbgEntry,String,String)> = Vec::new();
                for n2 in list_of_n2s {
                    //let n2_entry = dbg_nodes.get(&n2).unwrap(); // n2_entry didn't live long
                    //enough so i'm lamely copying the value. this is a TODO opt to pass it as
                    //reference
                    let n2_entry = dbg_nodes.get(n2).unwrap();
                    let mut vec_add_edge = | ori1 :&str,ori2 :&str|{
                        potential_edges.push((n2_entry.value().clone(),ori1.to_owned(),ori2.to_owned()));
                    };
                    let rev_n2 = n2.reverse();
                    if n1.suffix() == n2.prefix() {
                        vec_add_edge("+","+");
                    }
                    if revcomp_aware {
                        if n1.suffix() == rev_n2.prefix() {
                            vec_add_edge("+","-");
                        }
                        if rev_n1.suffix() == n2.prefix() {
                            vec_add_edge("-","+");
                        }
                        if rev_n1.suffix() == rev_n2.prefix() {
                            vec_add_edge("-","-");
                        }
                    }
                }
                if potential_edges.len() == 0 { continue;}
                let abundance_max = potential_edges.iter().map(|x| x.0.abundance).max().unwrap();
                let abundance_ref = std::cmp::min(abundance_max,n1_abundance);
                // write those edges to gfa
                for edge in potential_edges.iter()
                {
                    let (n2_entry,ori1,ori2) = edge;
                    let n2_abundance = n2_entry.abundance;
                    let n2_index     = n2_entry.index;
                    let n2_seqlen    = n2_entry.seqlen;
                    if presimp > 0.0 && potential_edges.len() >= 2
                    {
                        if (n2_abundance as f32) < presimp*(abundance_ref as f32) {
                            presimp_removed += 1;
                            removed_edges.insert((n1_index,n2_index));
                            continue;
                        }
                    }
                    let shift = if ori1 == "+" { n1_entry.shift.0 } else { n1_entry.shift.1 };
                    let overlap_length = std::cmp::min( n1_entry.seqlen - shift as u32, n2_entry.seqlen - 1 );
                    //if overlap_length == n2_entry.seqlen - 1 { println!("huh, node {} (minus shift {}) has overlap length as long as its neighbor {}",n1_seqlen,shift, n2_seqlen); }
                    if presimp == 0.0
                    { // no presimp means we can write edges now
                        let l_line = format!("L\t{}\t{}\t{}\t{}\t{}M\n", n1_index, ori1, n2_index, ori2, overlap_length);
                        write!(gfa_file, "{}", l_line).expect("error writing l_line");
                        nb_edges += 1;
                    }
                    else
                    { 
                        // write edges later, once we know if reverse wasn't affected by presimp
                        vec_edges.push((n1_index, ori1.clone(), n2_index, ori2.clone(), overlap_length));
                    }

                };
            }
        }
    }
    if presimp > 0.0
    { // write edges now because we couldn't earlier
        for edge in vec_edges.iter()
        {
            let (n1_index, ori1, n2_index, ori2, overlap_length) = edge;
            if removed_edges.contains(&(*n1_index,*n2_index)) ||removed_edges.contains(&(*n2_index,*n1_index))
            { // don't insert an edge if its reverse was deleted by presimp
                continue;
            }
            let l_line = format!("L\t{}\t{}\t{}\t{}\t{}M\n", *n1_index, *ori1, *n2_index, *ori2, overlap_length);
            write!(gfa_file, "{}", l_line).expect("error writing l_line");
            nb_edges += 1;
        }
    }
    println!("edges: {}", nb_edges);
    if presimp > 0.0 {
        println!("presimp={} removed {} edges",presimp,presimp_removed);
    }

    // create a real bidirected dbg object using petgraph
    /*{
      let mut gr = DiGraph::<Kmer,Kmer>::new();
      let mut node_indices : HashMap<Kmer,NodeIndex> = HashMap::new(); // bit redundant info, as nodes indices are in order of elements in dbg_nodes already; but maybe don't want to binary search inside it.
      for node in dbg_nodes.keys() { 
      let index = gr.add_node(node.clone());
      node_indices.insert(node.clone(),index);
      }
      let vec_edges : Vec<(NodeIndex,NodeIndex)> = dbg_edges.iter().map(|(n1,n2)| (node_indices.get(&n1).unwrap().clone(),node_indices.get(&n2).unwrap().clone())).collect();

      gr.extend_with_edges( vec_edges );
      }*/

    // graphml output
    //let graphml = GraphMl::new(&gr).pretty_print(true);
    //std::fs::write("graph.graphml", graphml.to_string()).unwrap();

    // graph pre-simplifications to avoid taking erroneous paths during naive coverage-oblivious
    // gfatools simplifications
    /*if presimp > 0.0
      {
      let removed_edges_all = presimp::find_removed_edges(&gr, &dbg_nodes, presimp, threads);
      presimp::presimp(&mut gr, &removed_edges_all);
      println!("{:?} edges removed during presimplification.", removed_edges_all.len());
      }*/ 

    // gfa output
    //println!("writing GFA..");
    //let mut node_indices = gfa_output::output_gfa(&gr, &dbg_nodes, &output_prefix, &kmer_seqs, &int_to_minimizer, levenshtein_minimizers, &node_indices);

    // write sequences of minimizers for each node
    // and also read sequences corresponding to those minimizers
    //println!("writing sequences..");
    //seq_output::write_minimizers_and_seq_of_kmers(&output_prefix, &mut node_indices, &kmer_origin, &dbg_nodes, k, l);
    //remove_file(seq_path);

    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Max RSS: {:?} GB", (get_memory_rusage() as f32)/1024.0/1024.0/1024.0);
}

