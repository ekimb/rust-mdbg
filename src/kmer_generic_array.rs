// notes:
// does not compile. at all. 
// i feel like i was close, but i don't know enough rust
// i couldn't define a type for arrays that are 1 element shorter

use generic_array::{ArrayLength,GenericArray};
use generic_array::sequence::Shorten;
use typenum::{Diff, U1};
use std::ops::Sub;
use arrayvec::ArrayVec;
use super::revcomp_aware;
use std::hash::{Hash, Hasher};

pub struct KmerArray<K> 
where K: ArrayLength<u32> + Sub<U1>
{
    data: GenericArray<u32, K>
}

impl <K> KmerArray<K>
where K: ArrayLength<u32> + Sub<U1>,
      Diff<K, U1>: ArrayLength<u32>
{

    pub fn suffix(&self) -> KmerArray<Diff<K,U1>>
    {
        let (init, last) = self.data.pop_back();
        init
    }
    
    pub fn prefix(&self) -> KmerArray<Diff<K,U1>>
    {
        // TODO adapt
        let (init, last) = self.data.pop_back();
        init
    }
    
    pub fn reverse_node(&self) -> KmerArray<K>
    {
        let rev_node_tmp :ArrayVec<KmerArray<K>> = &self.as_slice().into_iter().rev().map(|x| *x).collect();
        rev_node_tmp.into_inner().unwrap()
    }

    pub fn normalize(&self) -> KmerArray<K>
    {
        if !revcomp_aware { return self}
        std::cmp::min(*self,self.reverse())
    }

    pub fn make_from(ar: &[u32]) -> KmerArray<K>
    {
        KmerArray::<K> {
        //data: GenericArray::clone_from_slice(ar);
        }
        
    }



}

impl <K> PartialEq for KmerArray<K> 
where K: ArrayLength<u32>  + Sub<U1>
{
    fn eq(&self, other: &KmerArray<K>) -> bool {
        self.data == other.data
    }
}

impl <K> Eq for KmerArray<K> 
where K: ArrayLength<u32> + Sub<U1>
{
}

impl <K> Hash for KmerArray<K> 
where K: ArrayLength<u32> + Sub<U1>

{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.data.hash(state);
    }
}


