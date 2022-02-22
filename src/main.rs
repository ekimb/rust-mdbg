// rust-mdbg v0.1.0
// Copyright 2020-2021 Baris Ekim, Rayan Chikhi.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.
#[deny(clippy::mut_from_ref)]
//use petgraph::graph::NodeIndex;
//use petgraph_graphml::GraphMl;
//use petgraph::graph::DiGraph;
use pbr::ProgressBar;
use std::io::stderr;
//use std::error::Error;
use std::io::Write;
use std::io::{BufWriter, BufRead, BufReader};
use std::collections::HashMap;
use std::collections::hash_map::Entry;
use itertools::Itertools;
use closure::closure;
use crate::read::Read;
use std::collections::HashSet;
extern crate array_tool;
//use std::fs::remove_file;
use crossbeam_utils::{thread};
use structopt::StructOpt;
use std::sync::{Arc, Mutex};
use std::path::PathBuf;
use std::time::{Instant};
use std::fs;
use std::fs::File;
use std::mem::{MaybeUninit};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use lzzzz::lz4f::{WriteCompressor, BufReadDecompressor, Preferences};//, PreferencesBuilder, CLEVEL_HIGH};
use xx_bloomfilter::Bloom;
use flate2::read::GzDecoder;
use std::sync::atomic::{AtomicUsize, Ordering};
use glob::glob;
use dashmap::DashMap;
use std::cell::UnsafeCell;
use std::io::Result;
//use std::fmt::Arguments;
mod utils;
mod minimizers;
mod ec_reads;
mod kmer_vec;
mod poa;
mod read;
mod pairwise;

const REVCOMP_AWARE: bool = true; // shouldn't be set to false except for strand-directed data or for debugging
type Kmer = kmer_vec::KmerVec;
type Overlap = kmer_vec::KmerVec;
type DbgIndex = u32;// heavily optimized assuming we won't get more than 2B kminmers of abundance <= 65535
type DbgAbundance = u16;
type Repr = Vec<(Kmer, String, bool, String, (usize, usize))>;
type CorrMap = HashMap<String, (String, Vec<String>, Vec<usize>, Vec<u64>)>;

type Found = Option<(Repr, Read)>;
#[derive(Debug, Clone)] // seems necessary to move out of the Arc into dbg_nodes_view
struct DbgEntry {index: DbgIndex, abundance: DbgAbundance, seqlen: u32, shift: (u16, u16)} 
struct SeqFileType(WriteCompressor<File>);
unsafe impl Sync for SeqFileType {} // same trick as below. we won't share files among threads but Rust can't know that.
impl SeqFileType {
    fn new(v: File) -> SeqFileType {
        SeqFileType(WriteCompressor::new(v, Preferences::default()).unwrap())
    }
}
impl Write for SeqFileType {
    fn write(&mut self, buf: &[u8]) -> Result<usize> {
        self.0.write(buf)
    }

    fn flush(&mut self) -> Result<()> {
        self.0.flush()
    }
}
pub struct RacyBloom(UnsafeCell<Bloom>); // intentionnally allowing data races as a tradeoff for bloom speed
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
        unsafe {&mut *self.0.get()}
    }
}
type ThreadIdType = usize;
pub struct Params {
    l: usize,
    k: usize,
    n: usize,
    t: usize,
    density: f64,
    lmer_counts_min: u32,
    lmer_counts_max: u32,
    min_kmer_abundance: DbgAbundance,
    correction_threshold: i32,
    distance: usize,
    reference: bool,
    uhs: bool,
    lcp: bool,
    error_correct: bool,
    has_lmer_counts: bool,
    use_bf: bool,
    use_hpc: bool,
    use_syncmers: bool,
    s: usize,
    no_basespace: bool,
    debug: bool,
}


/*fn debug_output_read_minimizers(seq_str: &String, read_minimizers: &Vec<String>, read_minimizers_pos: &Vec<u32>) {
    println!("\nseq: {}", seq_str);
    print!("min: ");
    let mut current_minimizer :String = "".to_string();
    for i in 0..seq_str.len() {
        if read_minimizers_pos.contains(&(i as u32)) {
            let index = read_minimizers_pos.iter().position(|&r| r == i as u32).unwrap();
            current_minimizer = read_minimizers[index].clone();
            let c = current_minimizer.remove(0);
            if c == seq_str.chars().nth(i).unwrap() {print!("X");}
            else {print!("x");}
            continue;
        }
        if current_minimizer.len() > 0 {
            let c = current_minimizer.remove(0);
            print!("{}", c);
        }
        else {print!(".");}
    }
    println!("");
}*/

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
fn thread_update_hashmap<U, V>(hashmap_all: &Arc<Mutex<HashMap<usize, HashMap<U, V>>>>, hashmap: HashMap<U, V>, thread_num: usize) {
    let mut hashmap_all = hashmap_all.lock().unwrap();
    let entry = hashmap_all.entry(thread_num).or_insert_with(HashMap::new);
    *entry = hashmap; // I believe hashmap is moved in this function as per https://stackoverflow.com/a/29490907 
}

pub fn thread_update_vec<U>(vec_all: &Arc<Mutex<HashMap<usize, Vec<U>>>>, vec: Vec<U>, thread_num: usize) {
    let mut vec_all = vec_all.lock().unwrap();
    let entry = vec_all.entry(thread_num).or_insert_with(Vec::new);
    *entry = vec;
}

fn get_reader(path: &PathBuf) -> Box<dyn BufRead + Send> {
    let mut filetype = "unzip";
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
            Ok(file) => file,
            Err(error) => panic!("Error opening compressed file: {:?}.", error),
        };
    if filename_str.ends_with(".gz")  {filetype = "zip";}
    if filename_str.ends_with(".lz4") {filetype = "lz4";}
    let reader :Box<dyn BufRead + Send> = match filetype { 
        "zip" => Box::new(BufReader::new(GzDecoder::new(file))), 
        "lz4" => Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),
        _ =>     Box::new(BufReader::new(file)), 
    }; 
    reader
}

fn read_first_n_reads(filename: &PathBuf, fasta_reads: bool, max_reads: usize) -> (usize, usize) {
    let mut mean_length = 0;
    let mut max_length = 0;
    let mut nb_reads = 0;
    let buf = get_reader(&filename);
    // todo should factorize
    if fasta_reads {
        let mut reader = seq_io::fasta::Reader::new(buf);
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            let seq_str = String::from_utf8_lossy(record.seq()).to_string(); // might induce a copy? can probably be optimized (see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html)
            let l = seq_str.len();
            mean_length += l;
            max_length = std::cmp::max(l, max_length);
            nb_reads += 1;
            if nb_reads == max_reads {break;}
        }
    }
    else {
        let mut reader = seq_io::fastq::Reader::new(buf);
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            let seq_str = String::from_utf8_lossy(record.seq()).to_string(); // might induce a copy? can probably be optimized (see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html)
            let l = seq_str.len();
            mean_length += l;
            max_length = std::cmp::max(l, max_length);
            nb_reads += 1;
            if nb_reads == max_reads {break;}
        }
    }
    mean_length /= nb_reads;
    (mean_length, max_length) 
}

fn autodetect_k_l_d(filename: &PathBuf, fasta_reads: bool) -> (usize, usize, f64) {
    println!("Parsing input sequences to estimate mean read length...");
    let (mean_length, _max_length) = read_first_n_reads(&filename, fasta_reads, 100);
    println!("Detected mean read length of {} bp.",mean_length);
    // a bit crude, but let's try
    let d = 0.003;
    //let coeff : f64 = 3.0/4.0;
    let slightly_below_readlen : f64 = mean_length as f64;
    let k = (d * slightly_below_readlen) as usize;
    let l = 12;
    println!("Setting k = {} l = {} density = {}.", k, l, d);
    (k, l, d)
}

#[derive(Debug, StructOpt)]
#[structopt(name = "rust-mdbg",
    global_settings = &[structopt::clap::AppSettings::ColoredHelp, structopt::clap::AppSettings::ArgRequiredElseHelp])]
/// Original implementation of minimizer-space de Bruijn graphs (mdBG) for genome assembly.
///
/// rust-mdbg is an ultra-fast minimizer-space de Bruijn graph (mdBG) implementation, 
/// geared towards the assembly of long and accurate reads such as PacBio HiFi.
/// rust-mdbg is fast because it operates in minimizer-space, meaning that the reads, the assembly graph, 
/// and the final assembly, are all represented as ordered lists of minimizers, 
/// instead of strings of nucleotides. A conversion step then yields a classical base-space representation.
struct Opt {
    /// Activate debug mode
    ///
    /// Debug mode shows the base-space sequence and 
    /// the minimizer-space representation obtained from each read.
    #[structopt(long)]
    debug: bool,
    /// Input file (raw or gzip-/lz4-compressed FASTX)
    ///
    /// Input file can be FASTA/FASTQ, as well as gzip-compressed (.gz) or
    /// lz4-compressed (.lz4). Lowercase bases are currently not supported;
    /// see documentation for formatting.
    #[structopt(parse(from_os_str))]
    reads: PathBuf,
    #[structopt(long)]
    /// Universal k-mer file (enables universal hitting sets (UHS))
    ///
    /// Universal k-mers need to be provided as a single file,
    /// one universal k-mer per line. The minimizers will be selected
    /// if they are a universal k-mer and if they satisfy the density
    /// bound.
    uhs: Option<String>,
    #[structopt(long)]
    /// Core substring file (enables locally consistent parsing (LCP))
    ///
    /// Core substrings need to be provided as a single file,
    /// one core substring per line. The minimizers will be selected
    /// if they are a core substring and if they satisfy the density
    /// bound.
    lcp: Option<String>,
    /// Output prefix for GFA and .sequences files
    ///
    /// All files generated by rust-mdbg (including the GFA
    /// and .sequences output files) will have this file name.
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,
    /// k-min-mer length
    ///
    /// The length of each node of the mdBG. If
    /// fewer l-mers than this value are obtained
    /// from a read, they will be ignored.
    #[structopt(short, long)]
    k: Option<usize>,
    /// l-mer (minimizer) length
    ///
    /// The length of each minimizer selected using
    /// the minimizer scheme from base-space sequences.
    #[structopt(short, long)]
    l: Option<usize>,
    /// Tuple length for bucketing similar reads
    ///
    /// Reads that share a tuple (of l-mers) of this
    /// length will be bucketed as candidates for their
    /// respective partial order alignment (POA) graphs.
    #[structopt(short, long)]
    n: Option<usize>,
    /// POA path weight threshold
    ///
    /// During the consensus step for a partial order
    /// alignment (POA) graph, paths with weights below this value
    /// will be removed.
    #[structopt(short, long)]
    t: Option<usize>,
    /// Density threshold for density-based selection scheme
    /// 
    /// The density threshold is analogous to the
    /// fraction of l-mers that will be selected as
    /// minimizers from a read.
    #[structopt(short, long)]
    density: Option<f64>,
    /// Minimum k-min-mer abundance
    ///
    /// k-min-mers that occur fewer times than this
    /// value will be removed from the mdBG.
    #[structopt(long)]
    minabund: Option<usize>,
    /// Distance metric (0: Jaccard, 1: containment, 2: Mash)
    ///
    /// rust-mdbg uses a distance metric to filter out reads
    /// that are not very similar to a query read before constructing
    /// a partial order alignment (POA) graph.
    #[structopt(long)]
    distance: Option<usize>,
    /// POA correction threshold
    ///
    /// The maximum number of reads in a minimizer-space
    /// partial order alignment (POA) graph that can be
    /// replaced with the consensus.
    #[structopt(long)]
    correction_threshold: Option<i32>,
    /// Enable error correction with minimizer-space POA
    ///
    /// Partial order alignment (POA) can be used (in minimizer-space)
    /// to error-correct reads with up to 5% error rate.
    #[structopt(long)]
    error_correct: bool,
    /// Assemble error-corrected reads
    ///
    /// Reads that are error-corrected with rust-mdbg
    /// can be used as input.
    #[structopt(long)]
    restart_from_postcor: bool,
    /// Reference genome input
    ///
    /// Indicates that the input is a (single or a set of) 
    /// genome(s), not reads. Allows multi-line FASTA and
    /// doesn't filter any kminmers
    #[structopt(long)]
    reference: bool,
    /// Enable Bloom filters
    ///
    /// Bloom filters can be used to reduce memory usage,
    /// but results in slightly less contiguous assemblies.
    #[structopt(long)]
    bf: bool,
    /// Homopolymer-compressed (HPC) input
    ///
    /// Both raw and homopolymer-compressed (HPC) reads can
    /// be provided as input. If the reads are not compressed,
    /// rust-mdbg manually performs HPC, but uses the raw sequences
    /// for transformation into base-space.
    #[structopt(long)]
    hpc: bool,
    /// to save disk space, don't write the sequences in base-space
    /// corresponding to each k-min-mer
    #[structopt(long)]
    no_basespace: bool,
    /// use syncmers instead of universe minimizers
    #[structopt(long)]
    syncmers: bool,
    /// syncmer substring length
    #[structopt(short, long)]
    s: Option<usize>,
    /// l-mer counts (enables downweighting of frequent l-mers)
    /// use syncmers instead of universe minimizers
    ///
    /// Frequencies of l-mers in the reads (obtained using k-mer counters)
    /// can be provided in order to downweight frequently-occurring l-mers
    /// and increase contiguity.
    #[structopt(parse(from_os_str), long)]
    lmer_counts: Option<PathBuf>,
    /// Minimum l-mer count threshold
    ///
    /// l-mers with frequencies below this threshold will be
    /// downweighted.
    #[structopt(long)]
    lmer_counts_min: Option<u32>,
    /// Maximum l-mer count threshold
    ///
    /// l-mers with frequencies above this threshold will be
    /// downweighted.
    #[structopt(long)]
    lmer_counts_max: Option<u32>,
    /// Pre-simplification (pre-simp) threshold
    ///
    /// Additional graph simplification heuristics prior to the
    /// construction of the mdBG can be performed. The default
    /// value for the pre-simplification step is 0.1.
    #[structopt(long)]
    presimp: Option<f32>,
    /// Number of threads
    ///
    /// rust-mdbg is highly parallelized to decrease running
    /// time, but can be run on a single core as well.
    #[structopt(long)]
    threads: Option<usize>,
}

fn main() {
    let start = Instant::now();
    let opt = Opt::from_args();      
    let mut uhs : bool = false;
    let mut lcp : bool = false;
    let mut lmer_counts_filename = PathBuf::new();
    let mut uhs_filename = String::new();
    let mut lcp_filename = String::new();
    let mut output_prefix;
    let mut k : usize = 10;
    let mut l : usize = 12;
    let mut n : usize = 2;
    let mut t : usize = 0;
    let mut s : usize = 4; // syncmer mini-kmer size
    let mut density : f64 = 0.10;
    let mut min_kmer_abundance : DbgAbundance = 2;
    let mut distance : usize = 0;
    let mut error_correct : bool = false;
    let mut restart_from_postcor : bool = false;
    let mut correction_threshold : i32 = 0;
    let mut reference : bool = false;
    let mut has_lmer_counts : bool = false;
    let mut lmer_counts_min : u32 = 2;
    let mut lmer_counts_max : u32 = 100000;
    let mut presimp : f32 = 0.01;
    let mut use_bf : bool = false;
    let mut use_hpc : bool = false;
    let mut use_syncmers : bool = false;
    let mut no_basespace : bool = false;
    let mut threads : usize = 8;
    if opt.error_correct {error_correct = true;}
    if opt.reference {reference = true; error_correct = false;}
    let filename = opt.reads;
    if filename.as_os_str().is_empty() {panic!("Please specify an input file.");}

    let mut fasta_reads : bool = false;
    let filename_str = filename.to_str().unwrap();
        if filename_str.contains(".fasta.") || filename_str.contains(".fa.") || filename_str.ends_with(".fa") || filename_str.ends_with(".fasta") { // not so robust but will have to do for now
            fasta_reads = true;
            println!("Input file: {}", filename_str);
            println!("Format: FASTA");
    }
    if opt.k.is_none() && opt.l.is_none() && opt.density.is_none() {
        println!("Autodetecting values for k, l, and density.");
        let (ak, al, ad) = autodetect_k_l_d(&filename, fasta_reads);
        k = ak; l = al; density = ad;
    }
    else {
        if opt.k.is_some() {k = opt.k.unwrap()} else {println!("Warning: Using default k value ({}).", k);} 
        if opt.l.is_some() {l = opt.l.unwrap()} else {println!("Warning: Using default l value ({}).", l);}
        if opt.density.is_some() {density = opt.density.unwrap()} else {println!("Warning: Using default density value ({}%).", density * 100.0);}
    }
    if opt.n.is_some() {n = opt.n.unwrap()} else if error_correct {println!("Warning: Using default n value ({}).", n); }
    if opt.t.is_some() {t = opt.t.unwrap()} else if error_correct {println!("Warning: Using default t value ({}).", t); }
    if opt.minabund.is_some() {min_kmer_abundance = opt.minabund.unwrap() as DbgAbundance} else {println!("Warning: Using default minimum k-mer abundance value ({}).", min_kmer_abundance);}
    if opt.presimp.is_some() {presimp = opt.presimp.unwrap();} else {println!("Warning: Using default pre-simp value (0.01).");}
    if opt.threads.is_some() {threads = opt.threads.unwrap();} else {println!("Warning: Using default number of threads (8).");}
    if opt.correction_threshold.is_some() {correction_threshold = opt.correction_threshold.unwrap()} else if error_correct {println!("Warning: using default correction threshold value ({}).", correction_threshold);}
    if opt.distance.is_some() {distance = opt.distance.unwrap()}
    if distance > 2 {distance = 2;}
    let distance_type = match distance {0 => "jaccard", 1 => "containment", 2 => "mash", _ => "mash"};
    if opt.distance.is_none() && error_correct {println!("Warning: Using default distance metric ({}).", distance_type);}
    if opt.restart_from_postcor {restart_from_postcor = true;}
    if opt.bf {use_bf = true;}
    if opt.hpc {use_hpc = true;}
    if opt.syncmers {use_syncmers = true;}
    if use_syncmers
    {
        if opt.s.is_some() {s = opt.s.unwrap()} else {println!("Warning: Using default s value ({}).", s);}
    }
    if opt.no_basespace {no_basespace = true;}
    output_prefix = PathBuf::from(format!("graph-k{}-d{}-l{}", k, density, l));
    if opt.lmer_counts.is_some() { 
        has_lmer_counts = true;
        lmer_counts_filename = opt.lmer_counts.unwrap(); 
        if opt.lmer_counts_min.is_some() {lmer_counts_min = opt.lmer_counts_min.unwrap();} else {println!("Warning: Using default l-mer minimum count ({}).", lmer_counts_min);}
        if opt.lmer_counts_max.is_some() {lmer_counts_max = opt.lmer_counts_max.unwrap();} else {println!("Warning: Using default l-mer maximum count ({}).", lmer_counts_max);}
    } 
    if opt.uhs.is_some() { 
        uhs = true;
        uhs_filename = opt.uhs.unwrap(); 
    }
    if opt.lcp.is_some() { 
        lcp = true;
        lcp_filename = opt.lcp.unwrap(); 
    } 
    if opt.prefix.is_some() {output_prefix = opt.prefix.unwrap();} else {println!("Warning: Using default output prefix ({}).", output_prefix.to_str().unwrap());}
    let debug = opt.debug;
    let mut params = Params { 
        l,
        k,
        n,
        t,
        density,
        lmer_counts_min,
        lmer_counts_max,
        min_kmer_abundance,
        distance,
        correction_threshold,
        reference,
        uhs,
        lcp,
        error_correct,
        has_lmer_counts,
        use_bf,
        use_hpc,
        use_syncmers,
        s,
        no_basespace,
        debug,
    };
    // init some useful objects
    let mut nb_reads : u64 = 0;
    // get file size for progress bar
    let metadata = fs::metadata(&filename).expect("Error opening input file.");
    let file_size = metadata.len();
    let mut pb = ProgressBar::on(stderr(),file_size);
    let mut lmer_counts : HashMap<String, u32> = HashMap::new();

    if has_lmer_counts {
        let lmer_counts_file = match File::open(lmer_counts_filename) {
            Err(why) => panic!("Couldn't load l-mer counts: {}.", why.to_string()),
            Ok(lmer_counts_file) => lmer_counts_file,
        }; 
        let mut br = BufReader::new(lmer_counts_file);
        loop {
            let mut line = String::new();
            let new_line = |line: &mut String, br: &mut BufReader<File>| {line.clear(); br.read_line(line).ok();};
            if let Err(_e) = br.read_line(&mut line) {break;}
            if line.is_empty()                    {break;}
            let trimmed = line.trim().to_string();   
            let vec : Vec<String> = trimmed.split(' ').map(String::from).collect();
            let lmer = vec[0].to_string();
            let lmer_rev = utils::revcomp(&lmer);
            let lmer = if lmer > lmer_rev {lmer} else {lmer_rev}; //don't trust the kmer counter to normalize like we do
            let count = vec[1].parse::<u32>().unwrap();
            lmer_counts.insert(lmer, count);               
            new_line(&mut line, &mut br);
        }
    }
    let mut minimizer_to_int : HashMap<String,u64> = HashMap::new();
    let mut int_to_minimizer : HashMap<u64,String> = HashMap::new();
    // only need to initialize the minimizer_to_int / int_to_minimizer array if we do POA or use robust minimizers
    // they can be costly for k=14
    if has_lmer_counts || error_correct {
        let res = minimizers::minimizers_preparation(&mut params, &lmer_counts);
        minimizer_to_int = res.0;
        int_to_minimizer = res.1;
    }
    let (_mean_length, _max_length) = read_first_n_reads(&filename, fasta_reads, 10);

    // queue_len:
    // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
    // also: controls how many reads objects are buffered during fasta/fastq
    // parsing
    let queue_len = if _mean_length >= 1000000
                    {2} else {200}; // not too many seqs in cache if we're parsing reference genomes

    let mut uhs_bloom : RacyBloom = RacyBloom::new(Bloom::new(if use_bf {500_000_000} else {1}, 1_000_000_000_000_000));
    let mut lcp_bloom : RacyBloom = RacyBloom::new(Bloom::new(if use_bf {500_000_000} else {1}, 1_000_000_000_000_000));

    if params.uhs {uhs_bloom = minimizers::uhs_preparation(&mut params, &uhs_filename)}
    if params.lcp {lcp_bloom = minimizers::lcp_preparation(&mut params, &lcp_filename)}

    // dbg_nodes is a hash table containing (kmers -> (index,count))
    // it will keep only those with count > 1
    let dbg_nodes     : Arc<DashMap<Kmer, DbgEntry>> = Arc::new(DashMap::new()); // it's a Counter
    //let mut bloom : RacyBloom = RacyBloom::new(Bloom::new_with_rate(if use_bf {100_000_000} else {1}, 1e-7)); // a bf to avoid putting stuff into dbg_nodes too early
    let bloom         : RacyBloom = RacyBloom::new(Bloom::new(if use_bf {500_000_000} else {1}, 1_000_000_000_000_000)); // same bf but making sure we use only 1 hash function for speed
    static NODE_INDEX     : AtomicUsize = AtomicUsize::new(0); // associates a unique integer to each dbg node

    // correction stuff
    let mut buckets : HashMap<Vec<u64>, Vec<String>> = HashMap::new();
    //let corrected_map   = HashMap::<String, (String, Vec<String>, Vec<usize>, Vec<u64>)>::new(); // reduce runtime of POA by simultaneous correction of template and its aligned reads. keeps track of reads that have been already corrected
    let postcor_path        = PathBuf::from(format!("{}.postcor", output_prefix.to_str().unwrap()));
    let poa_path            = PathBuf::from(format!("{}.poa",     output_prefix.to_str().unwrap()));
    let mut reads_by_id     = HashMap::<String, Read>::new();

    // delete all previous sequence files
    for path in glob(&format!("{}*.sequences", output_prefix.to_str().unwrap())).expect("Failed to read glob pattern.") {
        let path = path.unwrap();
        let path = path.to_str().unwrap(); // rust really requires me to split the let statement in two..
        println!("Removing old sequences file: {}.", &path);
        let _res = fs::remove_file(path);
    }
    let seq_write = |file: &mut SeqFileType, s| {let _res = write!(file, "{}", s);};
    let sequences_files : Arc<DashMap<ThreadIdType, SeqFileType>> = Arc::new(DashMap::new());
    let create_sequences_file = |thread_id: ThreadIdType| -> SeqFileType {
        let seq_path = PathBuf::from(format!("{}.{}.sequences", output_prefix.to_str().unwrap(), thread_id));
        let file = match File::create(&seq_path) {
            Err(why) => panic!("Couldn't create file: {}.", why.to_string()),
            Ok(file) => file,
        };
        //let mut sequences_file = BufWriter::new(file);
        let mut sequences_file = SeqFileType::new(file); // regular lz4f
        //let mut sequences_file = WriteCompressor::new(&mut file, PreferencesBuilder::new().compression_level(CLEVEL_HIGH).build()).unwrap();  // too slow
        seq_write(&mut sequences_file, format!("# k = {}\n",k));
        seq_write(&mut sequences_file, format!("# l = {}\n",l));
        seq_write(&mut sequences_file, "# Structure of remaining of the file:\n".to_string());
        seq_write(&mut sequences_file, "# [node name]\t[list of minimizers]\t[sequence of node]\t[abundance]\t[origin]\t[shift]\n".to_string());
        sequences_file
    };

    let add_kminmer =|node: &Kmer, seq: Option<&str>, seq_reversed: &bool, origin: &str, shift: &(usize, usize), sequences_file: &mut SeqFileType, thread_id: usize, read_seq: Option<&str>, read_offsets: Option<(usize, usize, usize)>|
    {
        let mut previous_abundance; // for convenience, is the abundance of the kmer _before_ it was seen now
        let mut cur_node_index: DbgIndex = 0 as DbgIndex;
        let mut contains_key;

        // this code takes care of kminmers having abundances 1 and 2
        if use_bf && (!params.reference) && params.min_kmer_abundance > 1 {
            // this code discards abundances of 1 as we're not doing nested BFs (for now)
            //unsafe {
                // 1) check if the kminmer is in the BF
                if !bloom.get().check_and_add(&node)   { 
                    // 2) if not, it wasn't seen before so abundance has to be 1, and it gets inserted into BF
                    // for technical purposes, let's record the abundance before it's inserted (0) and not
                    // the true abundance (1)
                    //previous_abundance = 0;
                    return; // actually there is nothing more to do
                } else {
                    // 3) it was already in BF. meaning it was also inserted in dbg_nodes. we will get its true abundance later
                    previous_abundance = 1; 
                }
            //}
            contains_key = dbg_nodes.contains_key(node);
        }
        else { // old version without bf, just record abundance=1 kminmers in the hash table too
            contains_key = dbg_nodes.contains_key(node);
            if contains_key {
                previous_abundance = 1; // a sufficient placeholder. all we know is that the kminmer was previously seen so has abundance >= 1
            }
            else {
                cur_node_index = NODE_INDEX.fetch_add(1, Ordering::Relaxed) as DbgIndex;
                let lowprec_shift = (shift.0 as u16, shift.1 as u16);
                previous_abundance = 0;
                let seqlen = match seq {Some(read) => read.len() as u32, None => read_offsets.unwrap().2 as u32};
                // simulate the bf by inserting with abundance 0; will be incremented to 1 in
                // the code below
                dbg_nodes.insert(node.clone(), DbgEntry{index: cur_node_index, abundance: 0, seqlen, shift: lowprec_shift}); 
                contains_key = true;
            }
        }
        //println!("abundance: {}",abundance);
        if params.reference || previous_abundance >= 1 {
            // now record what we will save
            // or, record in the hash table anyway to save later
            let lowprec_shift = (shift.0 as u16, shift.1 as u16);
            if contains_key {
                let mut entry_mut = dbg_nodes.get_mut(node).unwrap();
                cur_node_index = entry_mut.index;
                previous_abundance = entry_mut.abundance;
                if previous_abundance == min_kmer_abundance - 1 {
                    let seqlen = match seq {Some(read) => read.len() as u32, None => read_offsets.unwrap().2 as u32};
                    entry_mut.seqlen = seqlen;
                    entry_mut.shift = lowprec_shift;
                }
                entry_mut.abundance += 1;
            }
            else { 
                cur_node_index = NODE_INDEX.fetch_add(1, Ordering::Relaxed) as DbgIndex;
                let seqlen = match seq {Some(read) => read.len() as u32, None => read_offsets.unwrap().2 as u32 };
                dbg_nodes.insert(node.clone(), DbgEntry{index: cur_node_index, abundance: previous_abundance+1, seqlen, shift: lowprec_shift}); 
            }
        }

        if params.reference || previous_abundance >= 1 || min_kmer_abundance == 1 {
            // now record what we will save
            if params.error_correct && thread_id != 0 {return;} // in error correct mode, only write sequences during second pass
            if previous_abundance == (min_kmer_abundance - 1) {
                // only save each kminmer once, exact at the time it passes the abundance filter.
                // (note that this doesnt enable to control which seq we save based on
                // median seq length, unfortunately)
                let seq = match seq {Some(read) => read, None => &read_seq.unwrap()[read_offsets.unwrap().0..read_offsets.unwrap().1]};
                let seq = if *seq_reversed {utils::revcomp(&seq)} else {seq.to_string()};
                let seq_line = format!("{}\t{}\t{}\t{}\t{}\t{:?}",cur_node_index, node.print_as_string(), seq, "*", origin, shift);
                if !params.no_basespace
                {
                    seq_write(sequences_file, format!("{}\n", seq_line));
                }
            }
        }
    };

    let add_kminmers = |vec: &Repr, thread_id: usize| {
        // determine to which sequence files to write
        if sequences_files.get(&thread_id).is_none() {
            sequences_files.insert(thread_id, create_sequences_file(thread_id));
        }

        for (node, seq, seq_reversed, origin, shift) in vec.iter() {
            let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
            add_kminmer(&node, Some(&seq), seq_reversed, &origin, shift, &mut sequences_file, thread_id, None, None);
        }
    };

    if !restart_from_postcor {
        let mut ec_file = BufWriter::new(File::create("/dev/null").unwrap());
        if error_correct || reference { // write this file even for an uncorrected reference, needed by evaluate_ec.py
            ec_file = ec_reads::new_file(&output_prefix); // reads before correction
        }
        
        // worker thread
        let process_read_aux = |seq_str: &[u8], seq_id: &str| -> Found {
            let thread_id :usize =  thread_id::get();
            let output : Repr = Vec::new();
            let seq = std::str::from_utf8(seq_str).unwrap(); // see https://docs.rs/seq_io/0.4.0-alpha.0/seq_io/fasta/index.html
            // those two next lines do a string copy of the read sequence, because in the case of a
            // reference, we'll need to remove newlines. also, that sequence will be moved to the
            // Read object. could be optimized later
            let seq_for_ref = if reference  {seq.replace("\n", "").replace("\r", "") // seq_io might return newlines in fasta seq 
                                            } else {String::new()};
            let seq = if reference {seq_for_ref} else {seq.to_string()}; 
            let read_obj = Read::extract(&seq_id, seq, &params, &minimizer_to_int, &uhs_bloom, &lcp_bloom);
            //println!("Received read in worker thread, transformed len {}", read_obj.transformed.len());


            // that's the non-optimized version
            /*if read_obj.transformed.len() > k {
                output = read_obj.read_to_kmers(&params);
            }
            add_kminmers(&output, thread_id); 
            */
            
            // determine to which sequence files to write
            if sequences_files.get(&thread_id).is_none() {
                sequences_files.insert(thread_id, create_sequences_file(thread_id));
            }
            let mut sequences_file = sequences_files.get_mut(&thread_id).unwrap();
            if read_obj.transformed.len() > k {
                // a copy of read_to_kmers except that we don't construct sequences until
                // absolutely necessarry
                for i in 0..(read_obj.transformed.len() - k + 1) {
                    let mut node : Kmer = Kmer::make_from(&read_obj.transformed[i..i+k]);
                    let mut seq_reversed = false;
                    if REVCOMP_AWARE { 
                        let (node_norm, reversed) = node.normalize(); 
                        node = node_norm;
                        seq_reversed = reversed;
                    } 
                    let origin = "*".to_string(); // uncomment the line below to track where the kmer is coming from (but not needed in production)
                    let minimizers_pos = &read_obj.minimizers_pos;
                    let position_of_second_minimizer = match seq_reversed {
                        true => minimizers_pos[i+k-1] - minimizers_pos[i+k-2],
                        false => minimizers_pos[i+1] - minimizers_pos[i]
                    };
                    let position_of_second_to_last_minimizer = match seq_reversed {
                        true => minimizers_pos[i+1] - minimizers_pos[i],
                        false => minimizers_pos[i+k-1] - minimizers_pos[i+k-2]
                    };
                    let shift = (position_of_second_minimizer, position_of_second_to_last_minimizer);
                    let read_offsets = (read_obj.minimizers_pos[i] as usize, (read_obj.minimizers_pos[i+k-1] as usize + l), (read_obj.minimizers_pos[i+k-1] + 1 - read_obj.minimizers_pos[i] + 1));
                    add_kminmer(&node, None, &seq_reversed, &origin, &shift, &mut sequences_file, thread_id, Some(&read_obj.seq), Some(read_offsets));
                }
            }

            if error_correct || reference {Some((output, read_obj))}
            else {None}
        };
        let process_read_fasta = |record: seq_io::fasta::RefRecord, found: &mut Found| {
            let seq_str = record.seq(); 
            let seq_id = record.id().unwrap().to_string();
            *found = process_read_aux(&seq_str, &seq_id);
        };
        let process_read_fastq = |record: seq_io::fastq::RefRecord, found: &mut Found| {
            let seq_str = record.seq(); 
            let seq_id = record.id().unwrap().to_string();
            *found = process_read_aux(&seq_str, &seq_id);
        };


        // parallel fasta parsing, with a main thread that writes to disk and populates hash tables
        let mut main_thread = |found: &Found| { // runs in main thread
            nb_reads += 1;
            //println!("Received read in main thread, nb kmers: {}", vec.len());
            let debug_only_display_read_and_minimizers = false;
            if debug_only_display_read_and_minimizers {
                // debug: just displays the read id and the list of minimizers
                let (_vec, read_obj) = found.as_ref().unwrap();
                println!("{} {}", &read_obj.id.to_string(), &read_obj.transformed.to_vec().iter().join(" ")); 
            }
            else if error_correct || reference {
                let (_vec, read_obj) = found.as_ref().unwrap();
                reads_by_id.insert(read_obj.id.to_string(), read_obj.clone());
                if read_obj.transformed.len() >= n {
                    ec_reads::record(&mut ec_file, &read_obj.id.to_string(), &read_obj.seq, &read_obj.transformed.to_vec(), &read_obj.minimizers, &read_obj.minimizers_pos);
                    for i in 0..read_obj.transformed.len()-n+1 {
                        let n_mer = utils::normalize_vec(&read_obj.transformed[i..i+n].to_vec());
                        let _entry = buckets.entry(n_mer).or_insert_with(Vec::<String>::new);
                        //entry.push(read_obj.id.to_string());
                    }
                }
            }
            // Some(value) will stop the reader, and the value will be returned.
            // In the case of never stopping, we need to give the compiler a hint about the
            // type parameter, thus the special 'turbofish' notation is needed.
            None::<()>
        };

        let buf = get_reader(&filename);
        println!("Parsing input sequences...");
        if fasta_reads {
            let reader = seq_io::fasta::Reader::new(buf);
            let _res = read_process_fasta_records(reader, threads as u32, queue_len, process_read_fasta, |_record, found| {main_thread(found)});
        }
        else {
            let reader = seq_io::fastq::Reader::new(buf);
            let _res = read_process_fastq_records(reader, threads as u32, queue_len, process_read_fastq, |_record, found| {main_thread(found)});
        }

        pb.finish_print("Converted reads to k-min-mers.");
        println!("Number of reads: {}", nb_reads);


        // this part will correct the reads, and dump them to disk
        if error_correct {
            ec_reads::flush(&mut ec_file);
            let chunks = ec_reads::load(&output_prefix);
            let mut chunk_length = 1;
            let ec_entries      = Arc::new(Mutex::new(HashMap::<usize, Vec<(String, String, Vec<u64>, Vec<String>, Vec<usize>)>>::new()));
            let poa_entries     = Arc::new(Mutex::new(HashMap::<usize, HashMap<String, Vec<String>>>::new()));
            let mut ec_file_poa     = ec_reads::new_file(&poa_path);      // POA debug info (which reads were recruited per template, I think. Baris can correct/confirm)
            let mut ec_file_postcor = ec_reads::new_file(&postcor_path);  // reads after correction
            if chunks.len() > threads {chunk_length = chunks.len()/threads+1;}
            thread::scope(|s| {
                let mut guards = Vec::with_capacity(threads);
                for (thread_num, chunk) in chunks.chunks(chunk_length).enumerate() {
                    let ec_entries = ec_entries.clone();
                    let poa_entries = poa_entries.clone();
                    let guard = s.spawn(closure!(move chunk, ref params, ref int_to_minimizer, ref buckets, ref reads_by_id, |_| {
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
                                /*read_obj.seq = corrected_map[&read_obj.id].0.to_string();
                                read_obj.minimizers = corrected_map[&read_obj.id].1.to_vec();
                                read_obj.minimizers_pos = corrected_map[&read_obj.id].2.to_vec();
                                read_obj.transformed = corrected_map[&read_obj.id].3.to_vec();*/
                            }
                            ec_entry.push((read_obj.id.to_string(), read_obj.seq.to_string(), read_obj.transformed.to_vec(), read_obj.minimizers.to_vec(), read_obj.minimizers_pos.to_vec()));
                        }
                        thread_update_vec(&ec_entries, ec_entry, thread_num);
                        thread_update_hashmap(&poa_entries, poa_map, thread_num);
                    }));
                    guards.push(guard);
                }
            }).unwrap();
            pb.finish_print("Done with correction.");
            let mut ec_entries = ec_entries.lock().unwrap();
            let mut poa_entries = poa_entries.lock().unwrap(); 
            for thread_num in 0..threads {
                let ec = ec_entries.entry(thread_num).or_insert_with(Vec::new);
                for tuple in ec.iter() {
                    ec_reads::record(&mut ec_file_postcor, &tuple.0, &tuple.1, &tuple.2, &tuple.3, &tuple.4);
                }
                let poa = poa_entries.entry(thread_num).or_insert_with(HashMap::new);
                for (temp, vec) in poa.iter_mut() {ec_reads::record_poa(&mut ec_file_poa, temp, vec.to_vec());}
            }
            ec_reads::flush(&mut ec_file_postcor);
            ec_reads::flush(&mut ec_file_poa);
        }
    }

    // now load the error corrected reads and construct the dbg from it
    // it's a little slow because read_to_kmers is called single threaded, there's room for
    // introducing multithreading here
    if error_correct || restart_from_postcor {
        dbg_nodes.clear(); // will be populated by add_kminmers()
        let chunks = ec_reads::load(&postcor_path);
        NODE_INDEX.store(0, Ordering::Relaxed);
        for ec_record in chunks.iter() {
            let mut read_obj = Read {id: ec_record.seq_id.to_string(), minimizers: ec_record.read_minimizers.to_vec(), minimizers_pos: ec_record.read_minimizers_pos.to_vec(), transformed: ec_record.read_transformed.to_vec(), seq: ec_record.seq_str.to_string(), corrected: false};
            if read_obj.transformed.len() > k { 
                let output = read_obj.read_to_kmers(&params);
                add_kminmers(&output, 0);
            }
        }
    }

    for mut sequences_file in sequences_files.iter_mut() {sequences_file.flush().unwrap();}

    // now DBG creation can start
    println!("Number of nodes before abundance filter: {}", dbg_nodes.len());
    dbg_nodes.retain(|_x, c| c.abundance >= (min_kmer_abundance as DbgAbundance));
    println!("Number of nodes after abundance filter: {}", dbg_nodes.len());
    let path = format!("{}{}", output_prefix.to_str().unwrap(),".gfa");
    let mut gfa_file = match File::create(&path) {
        Err(why) => panic!("Couldn't create {}: {}.", path, why.to_string()),
        Ok(file) => file,
    };
    writeln!(gfa_file, "H\tVN:Z:1.0").expect("Error writing GFA header.");
    // index k-1-mers
    let dbg_nodes_view = Arc::try_unwrap(dbg_nodes).unwrap().into_read_only();
    let mut km_index : HashMap<Overlap, Vec<&Kmer>> = HashMap::new(); 
    for (node, entry) in dbg_nodes_view.iter() {
        //let node = item.key();
        //let entry = item.value();
        // take this iteration opportunity to write S lines
        let length = entry.seqlen;
        let s_line = format!("S\t{}\t{}\tLN:i:{}\tKC:i:{}\n", entry.index, "*", length, entry.abundance);
        write!(gfa_file, "{}", s_line).expect("Error writing S line.");
        let first : Overlap = node.prefix().normalize().0;
        let second : Overlap = node.suffix().normalize().0;
        let mut insert_km = |key,val| {
            match km_index.entry(key) {
                Entry::Vacant(ent) => {ent.insert(vec![val]);},
                Entry::Occupied(mut ent) => {ent.get_mut().push(val);}
            }
        };
        insert_km(first, node);
        insert_km(second, node);
    }
    let mut nb_edges = 0;
    let mut presimp_removed = 0;
    let mut removed_edges : HashSet<(DbgIndex, DbgIndex)> = HashSet::new();
    let mut vec_edges = Vec::new(); // for presimp, keep a list of all edges to insert

    // create a vector of dbg edges (if we want to create a petgraph)
    // otherwise just output edges to gfa without keeping them in mem
    for (n1, n1_entry) in dbg_nodes_view.iter() {
        //let n1 = item.key();
        //let n1_entry = item.value();
        let rev_n1 = n1.reverse();
        let n1_abundance = n1_entry.abundance;
        let n1_index     = n1_entry.index;
        let n1_seqlen    = n1_entry.seqlen;

        // bit of a rust noob way to code this, because i'm not too familiar with types yet..
        let key1 = n1.suffix().normalize().0;
        let key2 = n1.prefix().normalize().0;
        for key in [key1, key2].iter() {
            if km_index.contains_key(&key) {
                let list_of_n2s : &Vec<&Kmer> = km_index.get(&key).unwrap();
                let mut potential_edges : Vec<(&DbgEntry, String, String)> = Vec::new();
                for n2 in list_of_n2s {
                    let n2_entry = dbg_nodes_view.get(n2).unwrap();
                    let mut vec_add_edge = |ori1: &str, ori2: &str|{
                        potential_edges.push((n2_entry, ori1.to_owned(), ori2.to_owned()));
                    };
                    let rev_n2 = n2.reverse();
                    if n1.suffix() == n2.prefix() {
                        vec_add_edge("+", "+");
                    }
                    if REVCOMP_AWARE {
                        if n1.suffix() == rev_n2.prefix() {
                            vec_add_edge("+", "-");
                        }
                        if rev_n1.suffix() == n2.prefix() {
                            vec_add_edge("-", "+");
                        }
                        if rev_n1.suffix() == rev_n2.prefix() {
                            vec_add_edge("-", "-");
                        }
                    }
                }
                if potential_edges.is_empty() {continue;}
                let abundance_max = potential_edges.iter().map(|x| x.0.abundance).max().unwrap();
                let abundance_ref = std::cmp::min(abundance_max, n1_abundance);
                // write those edges to gfa
                for edge in potential_edges.iter() {
                    let (n2_entry, ori1, ori2) = edge;
                    let n2_abundance = n2_entry.abundance;
                    let n2_index     = n2_entry.index;
                    let n2_seqlen    = n2_entry.seqlen;
                    if presimp > 0.0 && potential_edges.len() >= 2 && (n2_abundance as f32) < presimp * (abundance_ref as f32){
                        presimp_removed += 1;
                        removed_edges.insert((n1_index, n2_index));
                        continue;
                    }
                    let shift = if ori1 == "+" {n1_entry.shift.0} else {n1_entry.shift.1};
                    let overlap_length = std::cmp::min(n1_seqlen - shift as u32, n2_seqlen - 1);
                    //if overlap_length == n2_entry.seqlen - 1 { println!("huh, node {} (minus shift {}) has overlap length as long as its neighbor {}",n1_seqlen,shift, n2_seqlen); }
                    if presimp == 0.0 { // no presimp means we can write edges now
                        let l_line = format!("L\t{}\t{}\t{}\t{}\t{}M\n", n1_index, ori1, n2_index, ori2, overlap_length);
                        write!(gfa_file, "{}", l_line).expect("Error writing L line.");
                        nb_edges += 1;
                    }
                    else { 
                        // write edges later, once we know if reverse wasn't affected by presimp
                        vec_edges.push((n1_index, ori1.clone(), n2_index, ori2.clone(), overlap_length));
                    }
                };
            }
        }
    }
    if presimp > 0.0 { // write edges now because we couldn't earlier
        for edge in vec_edges.iter() {
            let (n1_index, ori1, n2_index, ori2, overlap_length) = edge;
            if removed_edges.contains(&(*n1_index, *n2_index)) || removed_edges.contains(&(*n2_index, *n1_index)) { // don't insert an edge if its reverse was deleted by presimp
                continue;
            }
            let l_line = format!("L\t{}\t{}\t{}\t{}\t{}M\n", *n1_index, *ori1, *n2_index, *ori2, overlap_length);
            write!(gfa_file, "{}", l_line).expect("Error writing L line.");
            nb_edges += 1;
        }
    }
    println!("Number of edges: {}", nb_edges);
    if presimp > 0.0 {
        println!("Pre-simp = {}: {} edges removed.",presimp,presimp_removed);
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

    // write sequences of minimizers for each node
    // and also read sequences corresponding to those minimizers
    //println!("writing sequences..");
    //seq_output::write_minimizers_and_seq_of_kmers(&output_prefix, &mut node_indices, &kmer_origin, &dbg_nodes, k, l);
    //remove_file(seq_path);

    let duration = start.elapsed();
    println!("Total execution time: {:?}", duration);
    println!("Maximum RSS: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);
}

