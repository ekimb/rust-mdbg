// Various helper functions used throughout.

pub fn revcomp(dna: &str) -> String {
	dna.chars()
        .rev()
        .map(switch_base)
        .collect::<String>()
}

fn switch_base(c: char) -> char {
    match c {
        'a' => 't',
        'c' => 'g',
        't' => 'a',
        'g' => 'c',
        'u' => 'a',
        'A' => 'T',
        'C' => 'G',
        'T' => 'A',
        'G' => 'C',
        'U' => 'A',
        _ => 'N'
    }
}

/// prints only the first 2 digits of each minimizer hash
pub fn pretty_minvec(seq: &[u64])  -> String {
    let mut s = String::new();
    for x in seq.iter() {
        s = format!("{}{} ", s, x.to_string()[..2].to_string());
    }
    s
}

/// for sequences that aren't kminmers (shorter or longer)
pub fn normalize_vec(seq: &[u64]) -> Vec<u64> {
    let mut seq_rev = seq.to_owned();
    seq_rev.reverse();
    std::cmp::min(seq_rev, seq.to_owned())
}
