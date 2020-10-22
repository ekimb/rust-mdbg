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
        s = format!("{}{} ",s,x.to_string()[..2].to_string());
    }
    s
}

/// for sequences that aren't kminmers (shorer or longer)
pub fn normalize_vec(seq :&Vec<u64>)  -> Vec<u64>
{
    let mut seq_rev = seq.clone();
    seq_rev.reverse();
    std::cmp::min(seq_rev, seq.clone())
    //seq.clone() // FIXME
}


pub fn median(numbers: &Vec<u32>) -> u32 {
    if numbers.len() == 1 {
        return numbers[0]
    }
    let mut numbers = numbers.clone();
    numbers.sort();
    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        //mean(&vec![numbers[mid - 1], numbers[mid]]) as i32
        numbers[mid-1] // we don't need to be that precise
    } else {
        numbers[mid]
    }
}
