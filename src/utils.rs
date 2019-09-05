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
