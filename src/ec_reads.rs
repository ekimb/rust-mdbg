// unfinished module

use std::fs::{File,remove_file};
use std::error::Error;
use std::io::Write;
use itertools::Itertools;
use std::io::{BufWriter, BufRead, BufReader};
use std::path::PathBuf;

pub fn make_filename(output_prefix: &PathBuf) -> String
{
    format!("{}.ec_data",output_prefix.to_str().unwrap())
}

pub fn new_file(output_prefix: &PathBuf) -> BufWriter<File> 
{
    let path = make_filename(output_prefix);
    let file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path, why.description()),
        Ok(file) => file,
    };
    let file = BufWriter::new(file);
    file
}

pub fn delete_file(output_prefix: &PathBuf)
{
    let path = make_filename(output_prefix);
    remove_file(path).unwrap();
}

pub fn record_poa(file: &mut BufWriter<File>, seq_id: &str, poa_ids: Vec<String>)
{
write!(file, "{}\t", seq_id).expect("error writing EC info");
write!(file, "{}\n", poa_ids.iter().join("\t")).expect("error writing EC info");


}
pub fn record(file: &mut BufWriter<File>, seq_id: &str, seq_str :&str, read_transformed: &Vec<u32>, read_minimizers: &Vec<String>, read_minimizers_pos: &Vec<u32>)
{
    write!(file, "{}\n", seq_id ).expect("error writing EC info");
    write!(file, "{}\n", seq_str).expect("error writing EC info");
    write!(file, "{}\n", read_transformed.iter().join(" ")).expect("error writing EC info");
    write!(file, "{}\n", read_minimizers.iter().join(" ")).expect("error writing EC info");
    write!(file, "{}\n", read_minimizers_pos.iter().join(" ")).expect("error writing EC info");
}
pub fn flush(file: &mut BufWriter<File>)
{
    file.flush().ok();
}

pub struct EcRecord
{
    pub seq_id  :String,
    pub seq_str :String,
    pub read_transformed: Vec<u32>,
    pub read_minimizers: Vec<String>,
    pub read_minimizers_pos: Vec<u32>
}

pub fn load(output_prefix: &PathBuf) -> Vec<EcRecord>
{
    let path = make_filename(output_prefix);
    let file = match File::open(path) {
        Err(why) => panic!("couldn't load ec file: {}", why.description()),
        Ok(file) => file,
    }; 
    let mut br = BufReader::new(file);
    let mut res : Vec<EcRecord> = vec![];
    let to_vec_int = |s: &String| -> Vec<u32> { s.trim().split(' ').map(|s| s.parse().unwrap()).collect::<Vec<u32>>() };
    let to_vec_str = |s: &String| -> Vec<String> { s.trim().split(' ').map(String::from).collect::<Vec<String>>() };

    loop
    {
        let mut line = String::new();
        let new_line = |line: &mut String, br :&mut BufReader<File>| { line.clear(); br.read_line(line).ok(); };
        if let Err(e) = br.read_line(&mut line) { break; }
        if line.len() == 0                      { break; }
        let seq_id  = line.trim().to_string();                    new_line(&mut line, &mut br);
        let seq_str = line.trim().to_string();                    new_line(&mut line, &mut br);
        let read_transformed: Vec<u32> = to_vec_int(&line);   new_line(&mut line, &mut br);
        let read_minimizers: Vec<String> = to_vec_str(&line); new_line(&mut line, &mut br);
        let read_minimizers_pos: Vec<u32> = to_vec_int(&line);
        res.push( EcRecord { seq_id, seq_str, read_transformed, read_minimizers, read_minimizers_pos } );
    }
    res
}
