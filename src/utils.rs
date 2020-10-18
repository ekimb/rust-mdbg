pub fn revcomp(dna: &str) -> String{
    // result vector
    let mut rdna: String = String::with_capacity(dna.len()); 

    // iterate through the input &str
    for c in dna.chars().rev() {
        rdna.push(switch_base(c))
    }
    rdna
}

fn switch_base(c:char) -> char {
    match c {
        'a' => 't' ,
        'c' => 'g' ,
        't' => 'a' ,
        'g' => 'c' ,
        'u' => 'a',
        'A' => 'T' ,
        'C' => 'G' ,
        'T' => 'A' ,
        'G' => 'C',
        'U' => 'A',
        _ => 'N'
    }
}

/// prints only the first 2 digits of each minimizer hash
pub fn pretty_minvec(seq :&Vec<u64>)  -> String
{
    let mut s = String::new();
    for x in seq.iter()
    {
        s = format!("{}{} ",s,x.to_string()[..].to_string());
    }
    s
}
