// read_stats.rs
// Contains the "ReadStats" struct, to collect statistics of kminmer count for each read

use std::{fs::File, io::Write, mem::MaybeUninit, sync::Mutex};

static ENABLED: bool = true;

// some rust black magic to have multithreaded static file handle
// from https://stackoverflow.com/a/72187853
static mut STATS_FILE: MaybeUninit<Mutex<File>> = MaybeUninit::uninit();

#[derive(Clone, Debug, PartialEq)]
pub struct ReadStats {
    pub q_id: String,
    pub counts: Vec<u32>,
}

impl ReadStats {

    pub fn init(output_prefix: &str)
    {
        if ENABLED
        {
            // stats file generation
            let stats_path = format!("{}{}", output_prefix, ".read_stats");
            unsafe
            {
            STATS_FILE = match File::create(&stats_path) {
                Err(_why) => panic!("Couldn't create {}", stats_path),
                Ok(stats_file) =>  MaybeUninit::new(Mutex::new(stats_file)),
            };
            }
            println!("Stats module initialized.");
        }
    }

    // A new Stats object for a read
    pub fn new(q_id: &str) -> Self {
        let s = ReadStats {
            q_id: q_id.to_string(),
            counts: vec![]
        };
        s
    }
    
    // Add a new potential reference location
    pub fn add(&mut self, count: u32)
    {
        self.counts.push(count);
    }

    // Compute number of jumps and write to stats file
    pub fn finalize(&mut self)
    {
        let mut counts_str = String::new();
        for c in &self.counts {
            counts_str.push_str(&format!("{} ",c));
        }
        let stats_line = format!("{}: {}",self.q_id,counts_str);
        unsafe
        {
        write!(STATS_FILE.assume_init_mut().lock().unwrap(), "{}\n", stats_line).expect("Error writing stats line.");
        }
    }

}
