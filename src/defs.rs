use super::Params;
use super::minimizers;
use std::collections::HashMap;
#[derive(Clone, Debug)]
pub struct Read {
    pub seq: String,
    pub mins: Vec<u64>
}
impl PartialEq for Read
{
    fn eq(&self, other: &Read) -> bool {
        self.mins == other.mins
    }
}
impl  Read {
    pub fn new(seq : String, params : &Params) -> Self {
        let mins = minimizers::minhash(seq.as_bytes(), &params);
        Read {seq, mins}
    }
    
}
pub struct Bucket {
    ID: Vec<u64>,
    contents: Vec<Read>,
}
impl Bucket {
    pub fn new(ID: Vec<u64>) -> Self {
        Bucket {ID, contents : Vec::<Read>::new()}
    }
    pub fn insert(&mut self, read : Read) {
        self.contents.push(read);
    }
    pub fn contains(&self, read : Read) -> bool {
        self.contents.contains(&read)
    }
}
pub struct BucketCollection {
    pub map: HashMap<Vec<u64>, Bucket>,
}
impl BucketCollection {
    pub fn new() -> Self {
        BucketCollection {map : HashMap::new()}
    }
    pub fn get_contents(&mut self, ID: Vec<u64>) -> &Vec<Read>{
        &self.map[&ID].contents
    }
    pub fn exists(&mut self, ID: Vec<u64>) -> bool {
        self.map.contains_key(&ID)
    }
    pub fn populate_with_read(&mut self, read: Read, params: &Params) {
        let n = params.n;
        
    }

}
