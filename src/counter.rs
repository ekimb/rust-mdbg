//https://github.com/sstadick/kmer-cnt/blob/master/kc/src/main.rs

use fastq::{parse_path_fa, Record};
use fnv::FnvHashMap;
use std::env::args;
use std::io::{self};

type Counter = FnvHashMap<u64, usize>;

static SEQ_NT4_TABLE: &'static [u8] = &[
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
];

fn count_kmer(counter: &mut Counter, k: usize, seq: &[u8]) {
    let mut strands = (0, 0);
    let mask = (1u64 << k * 2) - 1;
    let shift = (k - 1) * 2;
    let mut i = 0;
    let mut l = 0;
    loop {
        if !(i < seq.len()) {
            break;
        }
        let c = SEQ_NT4_TABLE[seq[i] as usize];
        if c < 4 {
            // Not an N base
            strands.0 = (strands.0 << 2 | c as u64) & mask; // forward strand
            strands.1 = strands.1 >> 2 | (3 - c as u64) << shift; // reverse strand
            l += 1;
            if l >= k {
                // we find a kmer
                let y = if strands.0 < strands.1 {
                    strands.0
                } else {
                    strands.1
                };
                // only add one strand! TODO: why?
                let count = counter.entry(y).or_insert(0);
                *count += 1;
            }
        } else {
            // if there is an N, restart
            l = 0;
            strands.0 = 0;
            strands.1 = 0;
        }

        i += 1;
    }
}

fn print_hist(counter: &Counter) {
    let mut hist = vec![0; 256];
    for (_kmer, count) in counter {
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

fn main() -> io::Result<()> {
    // Get the file from
    let mut counter = Counter::default();
    let filename = args().nth(1);
    let path = match filename.as_ref().map(String::as_ref) {
        None | Some("-") => None,
        Some(name) => Some(name),
    };
    parse_path_fa(path, |parser| {
        parser
            .each(|rec| {
                count_kmer(&mut counter, 31, rec.seq());
                true
            })
            .unwrap();
    })?;
    print_hist(&counter);
    std::mem::forget(counter);
    Ok(())
}