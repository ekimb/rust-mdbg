// Copyright 2017-2018 Brett Bowman, Jeff Knaggs
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Partial-Order Alignment for fast alignment and consensus of multiple homologous sequences.
//!
//! For the original concept and theory, see:
//! * Lee, Christopher, Catherine Grasso, and Mark F. Sharlow. "Multiple sequence alignment using
//! partial order graphs." Bioinformatics 18.3 (2002): 452-464.
//! * Lee, Christopher. "Generating consensus sequences from partial order multiple sequence
//! alignment graphs." Bioinformatics 19.8 (2003): 999-1008.
//!
//! For a modern reference implementation, see poapy:
//! https://github.com/ljdursi/poapy
//!

use std::cmp::{max, Ordering};
use std::collections::HashMap;
use super::Params;
use petgraph::graph::NodeIndex;
//use petgraph::graph::EdgeIndex;
use petgraph::visit::EdgeRef;
use petgraph::algo::toposort;
use petgraph::visit::Topo;
//use petgraph::visit::Dfs;
use petgraph::dot::{Dot};//, Config};
use petgraph::{Directed, Graph, Incoming, Outgoing};
// for Smith-Waterman
use super::pairwise;
use super::utils::pretty_minvec;
pub type POAGraph = Graph<u64, (i32, String), Directed, usize>;

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.

pub const MIN_SCORE: i32 = -858_993_459;

/// Trait required to instantiate a Scoring instance
pub trait MatchFunc {
    fn score(&mut self, a: u64, b: u64) -> i32;
}

/// A concrete data structure which implements trait MatchFunc with constant
/// match and mismatch scores
#[derive(Debug, Clone)]
pub struct MatchParams {
    pub match_score: i32,
    pub mismatch_score: i32,
}

impl MatchParams {
    /*/// Create new MatchParams instance with given match and mismatch scores
    ///
    /// # Arguments
    ///
    /// * `match_score` - the score for a match (should not be negative)
    /// * `mismatch_score` - the score for a mismatch (should not be positive)
    pub fn new(match_score: i32, mismatch_score: i32) -> Self {
        assert!(match_score >= 0, "match_score can't be negative");
        assert!(mismatch_score <= 0, "mismatch_score can't be positive");
        MatchParams {
            match_score,
            mismatch_score,
        }
    }*/
}

impl MatchFunc for MatchParams {
    #[inline]
    fn score(&mut self, a: u64, b: u64) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The trait Matchfunc is also implemented for Fn(u8, u8) -> i32 so that Scoring
/// can be instantiated using closures and custom user defined functions
impl<F> MatchFunc for F
where
    F: FnMut(u64, u64) -> i32,
{
    fn score(&mut self, a: u64, b: u64) -> i32 {
        (self)(a, b)
    }
}

/// Details of scoring are encapsulated in this structure.
/// An affine gap score model is used so that the gap score for a length 'k' is:
/// GapScore(k) = gap_open + gap_extend * k
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct Scoring<F: MatchFunc> {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub match_fn: F,
    pub match_scores: Option<(i32, i32)>,
    pub xclip_prefix: i32,
    pub xclip_suffix: i32,
    pub yclip_prefix: i32,
    pub yclip_suffix: i32,
}

impl Scoring<MatchParams> {
    /*/// Create new Scoring instance with given gap open, gap extend penalties
    /// match and mismatch scores. The clip penalties are set to MIN_SCORE by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_score` - the score for a match
    /// * `mismatch_score` - the score for a mismatch
    pub fn from_scores(gap_open: i32, gap_extend: i32, match_score: i32, mismatch_score: i32) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn: MatchParams::new(match_score, mismatch_score),
            match_scores: Some((match_score, mismatch_score)),
            xclip_prefix: MIN_SCORE,
            xclip_suffix: MIN_SCORE,
            yclip_prefix: MIN_SCORE,
            yclip_suffix: MIN_SCORE,
        }
    }*/
}

impl<F: MatchFunc> Scoring<F> {
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// and the score function. The clip penalties are set to MIN_SCORE by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions (also see bio::scores)
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Scoring {
            gap_open,
            gap_extend,
            match_fn,
            match_scores: None,
            xclip_prefix: MIN_SCORE,
            xclip_suffix: MIN_SCORE,
            yclip_prefix: MIN_SCORE,
            yclip_suffix: MIN_SCORE,
        }
    }

    /*/// Sets the prefix and suffix clipping penalties for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Clipping penalty for x (both prefix and suffix, should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip(-5);
    /// assert!(scoring.xclip_prefix == -5);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == -5);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn xclip(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_prefix = penalty;
        self.xclip_suffix = penalty;
        self
    }

    /// Sets the prefix clipping penalty for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Prefix clipping penalty for x (should not be positive)
    ///
    /// # Example
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip_prefix(-5);
    /// assert!(scoring.xclip_prefix == -5);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn xclip_prefix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_prefix = penalty;
        self
    }

    /// Sets the suffix clipping penalty for x to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Suffix clipping penalty for x (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).xclip_suffix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == -5);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn xclip_suffix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.xclip_suffix = penalty;
        self
    }

    /// Sets the prefix and suffix clipping penalties for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Clipping penalty for y (both prefix and suffix, should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == -5);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == -5);
    /// ```
    pub fn yclip(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self.yclip_suffix = penalty;
        self
    }

    /// Sets the prefix clipping penalty for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Prefix clipping penalty for y (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip_prefix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == -5);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == MIN_SCORE);
    /// ```
    pub fn yclip_prefix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self
    }

    /// Sets the suffix clipping penalty for y to the input value
    ///
    /// # Arguments
    ///
    /// * `penalty` - Suffix clipping penalty for y (should not be positive)
    ///
    /// ```rust
    /// use bio::alignment::pairwise::{Scoring, MIN_SCORE};
    /// let scoring = Scoring::from_scores(0, -2, 1, -2).yclip_suffix(-5);
    /// assert!(scoring.xclip_prefix == MIN_SCORE);
    /// assert!(scoring.yclip_prefix == MIN_SCORE);
    /// assert!(scoring.xclip_suffix == MIN_SCORE);
    /// assert!(scoring.yclip_suffix == -5);
    /// ```
    pub fn yclip_suffix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_suffix = penalty;
        self
    }*/
}

#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}

pub struct Alignment {
    pub score: i32,
    ystart: usize,
    //    xstart: Edge,
    operations: Vec<AlignmentOperation>,
}

#[derive(Debug, Clone)]
pub struct TracebackCell {
    score: i32,
    op: AlignmentOperation,
}

impl Ord for TracebackCell {
    fn cmp(&self, other: &TracebackCell) -> Ordering {
        self.score.cmp(&other.score)
    }
}

impl PartialOrd for TracebackCell {
    fn partial_cmp(&self, other: &TracebackCell) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for TracebackCell {
    fn eq(&self, other: &TracebackCell) -> bool {
        self.score == other.score
    }
}

//impl Default for TracebackCell { }

impl Eq for TracebackCell {}

pub struct Traceback {
    rows: usize,
    cols: usize,

    // store the last visited node in topological order so that
    // we can index into the end of the alignment when we backtrack
    last: NodeIndex<usize>,
    matrix: Vec<Vec<TracebackCell>>,
}

impl Traceback {
    /// Create a Traceback matrix with given maximum sizes
    ///
    /// # Arguments
    ///
    /// * `m` - the number of nodes in the DAG
    /// * `n` - the length of the query sequence
    ///
    fn with_capacity(m: usize, n: usize) -> Self {
        let matrix = vec![
            vec![
                TracebackCell {
                    score: 0,
                    op: AlignmentOperation::Match(None)
                };
                n + 1
            ];
            m + 1
        ];
        Traceback {
            rows: m,
            cols: n,
            last: NodeIndex::new(0),
            matrix,
        }
    }

    /*/// Populate the edges of the traceback matrix
    fn initialize_scores(&mut self, gap_open: i32, gap_extend: i32) {
        for (i, row) in self
            .matrix
            .iter_mut()
            .enumerate()
            .take(self.rows + 1)
            .skip(1)
        {
            // TODO: these should be -1 * distance from head node
            row[0] = TracebackCell {
                score: gap_open + ((i-1) as i32) * gap_extend, 
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..=self.cols {
            self.matrix[0][j] = TracebackCell {
                score: gap_open + ((j-1) as i32) * gap_extend,
                op: AlignmentOperation::Ins(None),
            };
        }
    }*/

    fn new() -> Self {
        Traceback {
            rows: 0,
            cols: 0,
            last: NodeIndex::new(0),
            matrix: Vec::new(),
        }
    }

    fn set(&mut self, i: usize, j: usize, cell: TracebackCell) {
        self.matrix[i][j] = cell;
    }

    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        &self.matrix[i][j]
    }

    /*pub fn print(&self, g: &POAGraph, query: Vec<u64>) {
        let (m, n) = (g.node_count(), query.len());
        print!(".\t");
        for base in query.iter().take(n) {
            print!("{:?}\t", *base);
        }
        for i in 0..m {
            print!("\n{:?}\t", g.raw_nodes()[i].weight);
            for j in 0..n {
                print!("{}.\t", self.get(i + 1, j + 1).score);
            }
        }
        println!();
    }*/
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
pub struct Aligner<F: MatchFunc> {
    pub traceback: Traceback,
    query: Vec<u64>,
    query_str: String, 
    query_minim_pos: Vec<usize>,
    pub poa: Poa<F>,
}

impl<F: MatchFunc> Aligner<F> {
    /// Create new instance.
    pub fn new(scoring: Scoring<F>, reference: &[u64], ref_str: Option<&String>, ref_minim_pos: Option<&Vec<usize>>) -> Self {
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            query_str: ref_str.unwrap().clone(),
            query_minim_pos: ref_minim_pos.unwrap().clone(),
            poa: Poa::from_string(scoring, reference, ref_str, ref_minim_pos),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self) -> &mut Self {
        let alignment = self.alignment();
        self.poa.add_alignment(&alignment, self.query.to_vec(), Some(&self.query_str), Some(&self.query_minim_pos));
        self
    }

    /// Return alignment of last added query against the graph.
    pub fn alignment(&self) -> Alignment {
        let debug = false;

        let mut ops: Vec<AlignmentOperation> = vec![];
    
        // Now backtrack through the matrix to construct an optimal path
        let mut i;
        let mut j;

        // semiglobal
        let cands: Vec<NodeIndex<usize>> = self.poa.graph.node_indices().filter(|x| self.poa.graph.neighbors_directed(*x, Outgoing).count() == 0).collect();
        //println!("Cands {:?}, len {}", cands, cands.len());
        let max_cand = cands.iter().map(|cand| ((cand.index() + 1, self.traceback.cols), self.traceback.matrix[cand.index()+1][self.traceback.cols].score)).max_by(|a, b| a.1.partial_cmp(&b.1).unwrap()).unwrap();
        let i_0 = (max_cand.0).0;
        let j_0 = (max_cand.0).1;

        i = i_0;
        j = j_0;
        if debug{print!("traceback: ");}
        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.traceback.matrix[i][j].op.clone());
            match self.traceback.matrix[i][j].op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                    if debug{print!("M");}
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                    if debug{print!("-");}
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j -= 1; 
                    if debug{print!("I");}
                }
                AlignmentOperation::Match(None) => {
                    if debug{print!("m");}
                    j -= 1; //?
                    break;
                }
                AlignmentOperation::Del(None) => {
                    if debug{print!("d");}
                    break; //semiglobal can end anywhere in the first column
                    //j -= 1;
                }
                AlignmentOperation::Ins(None) => {
                    if debug{print!("i");}
                    i -= 1;
                    j -= 1; //?
                }
            }
        }
        if debug{println!();}

        ops.reverse();

        if debug{println!("ystart: {}",j);}
        Alignment {
            score: self.traceback.matrix[i_0][j_0].score,
            ystart: j,
            operations: ops,
        }
    }

    /*/// Globally align a given query against the graph.
    pub fn global(&mut self, query: &Vec<u64>) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.global(&query);
        self
    }*/

    /// Semiglobally align a given query against the graph.
    pub fn semiglobal(&mut self, query: &[u64], query_str: Option<&String>, query_minim_pos: Option<&Vec<usize>>) -> &mut Self {
        self.query = query.to_vec();
        self.query_str = query_str.unwrap().clone();
        self.query_minim_pos = query_minim_pos.unwrap().clone();
        self.traceback = self.poa.semiglobal(&query);
        self
    }

    /*/// Return alignment graph.
    pub fn graph(&self) -> &POAGraph {
        &self.poa.graph
    }*/

    /// print a pretty alignment of last query to the graph
    pub fn print_aln(&mut self) ->String {
        let alignment = self.alignment();
        self.poa.pretty(&alignment, &self.query.to_vec())
    }

    /// gets template boundary in consensus
    /// Perform Smith-Waterman on the consensus with the original template
    /// to determine the boundary of the template on the consensus
    pub fn consensus_boundary(&mut self, cns : &[u64], cns_es: &[String], orig: &[u64], debug: bool)
        -> (Vec<u64>, Vec<String>)
    {
        let x = &orig[..];
        let y = & cns[..];
        let score = |a: u64, b: u64| if a == b {1i32} else {-1i32};
        let mut regular_aligner = pairwise::Aligner::with_capacity(x.len(), y.len(), -1, -1, &score);
        let alignment = regular_aligner.semiglobal(x, y);
        if alignment.yend - alignment.ystart < 2 { return (Vec::new(), Vec::new()); }
        let adjusted_cns = cns   [alignment.ystart..alignment.yend].to_vec();
        let adjusted_es  = cns_es[alignment.ystart..(alignment.yend-1)].to_vec();
        if debug
        {
            println!("alignment: {:?} (score: {})",alignment.operations, alignment.score);
            println!("consensus    : {:?}",pretty_minvec(cns));
            println!("orig template: {:?}",pretty_minvec(orig));
            println!("adjusted cns : {:?}",pretty_minvec(&adjusted_cns));
        }

        // TODO comment that part if the warning below is never printed

        let mut y_rev = cns.to_owned();
        y_rev.reverse();
        let rev_alignment = regular_aligner.semiglobal(x, &y_rev);
        if false // debug
        {
            println!("bwd alignment: {:?} (score: {})",rev_alignment.operations, alignment.score);
        }
        if rev_alignment.score > alignment.score
        {
            println!("WARNING! (rev alignment POA consensus score {} higher than fwd alignment {})", rev_alignment.score, alignment.score);
        }

        (adjusted_cns, adjusted_es)
    }
 

}

/// A partially ordered alignment graph
///
/// A directed acyclic graph datastructure that represents the topology of a
/// traceback matrix.
///
pub struct Poa<F: MatchFunc> {
    scoring: Scoring<F>,
    pub graph: POAGraph,
}

impl<F: MatchFunc> Poa<F> {
    /*/// Create a new aligner instance from the directed acyclic graph of another.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `poa` - the partially ordered reference alignment
    ///
    pub fn new(scoring: Scoring<F>, graph: POAGraph) -> Self {
        Poa { scoring, graph}
    }*/

    /// Create a new POA graph from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    ///
    pub fn from_string(scoring: Scoring<F>, seq: &[u64], seq_str: Option<&String>, seq_minim_pos: Option<&Vec<usize>>) -> Self {
        let mut graph: POAGraph =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        let mut node: NodeIndex<usize>;
        let mut i = 1;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);

            // adding sequence between minimizers to graph edges
            let seq_minim_pos = seq_minim_pos.unwrap();
            let between_minims = &seq_str.unwrap()[seq_minim_pos[i-1] as usize..seq_minim_pos[i] as usize];

            graph.add_edge(prev, node, (1,String::from(between_minims)));
            //let edge = graph.find_edge(prev, node).unwrap();
            prev = node;
            i += 1;
        }

        Poa { scoring, graph }
    }
   
    /// determine whether the gap penalty is an open or an extend based on previous cell
    fn determine_gap_penalty(&self, traceback_cell : &TracebackCell, current_op: AlignmentOperation) -> i32
    {
        match traceback_cell.op {
            AlignmentOperation::Match(None) => {
                self.scoring.gap_open
            }
            AlignmentOperation::Match(Some((_, _p))) => {
                self.scoring.gap_open
            }
            AlignmentOperation::Ins(None) => {
                match current_op
                {
                    AlignmentOperation::Ins(None) =>  {
                        self.scoring.gap_extend
                    },
                     AlignmentOperation::Ins(Some(_)) => {
                        self.scoring.gap_extend
                    },
                    _ => {
                        self.scoring.gap_open
                    }
                }
            }
            AlignmentOperation::Ins(Some(_)) => {
                match current_op
                {
                    AlignmentOperation::Ins(None) =>  {
                        self.scoring.gap_extend
                    },
                     AlignmentOperation::Ins(Some(_)) => {
                        self.scoring.gap_extend
                    },
                    _ => {
                        self.scoring.gap_open
                    }
                }
            }
            AlignmentOperation::Del(_) => {
                match current_op
                {
                    AlignmentOperation::Del(None) =>  {
                        self.scoring.gap_extend
                    },
                    _ => {
                        self.scoring.gap_open
                    }
                }
            }
        }
    }

    /*/// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    ///

    pub fn global(&mut self, query: &Vec<u64>) -> Traceback {
        assert!(self.graph.node_count() != 0);

        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.scoring.gap_open, self.scoring.gap_extend);

        traceback.set(
            0,
            0,
            TracebackCell {
                score: 0,
                op: AlignmentOperation::Match(None),
            },
        );

        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            // query base and its index in the DAG (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    TracebackCell {
                        score: traceback.get(0, j - 1).score + self.scoring.match_fn.score(r, *q),
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        score: MIN_SCORE,
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        let gap_penalty = self.determine_gap_penalty(traceback.get(i_p, j), AlignmentOperation::Del(Some((i_p - 1, i))));
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + self.scoring.match_fn.score(r, *q),
                                        op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + gap_penalty,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };

                let gap_penalty = self.determine_gap_penalty(traceback.get(i, j-1),AlignmentOperation::Ins(Some(i - 1)));
                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + gap_penalty,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
            }
        }

        traceback
    }*/

    // code not ready yet. initialization is okay, but the score update + the traceback stopping condition (not in this function) aren't changed
     
    ///  semi-global alignment 
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    ///

    pub fn semiglobal(&mut self, query: &[u64]) -> Traceback {
        assert!(self.graph.node_count() != 0);

        // dimensions of the traceback matrix
        let (length, width) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(length, width);
        //traceback.initialize_scores(self.scoring.gap_open);
        for (_i, row) in traceback 
            .matrix
            .iter_mut()
            .enumerate()
            .take(traceback.rows + 1)
            .skip(1)
        {
            row[0] = TracebackCell {
                score: 0, //semi-global, we may start anywhere in the graph
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..=traceback.cols {
            traceback.matrix[0][j] = TracebackCell {
                score: (j as i32) * self.scoring.gap_open, // but if we consume j chars of the query, we pay j gap penalties
                op: AlignmentOperation::Ins(None),
            };
        }

        traceback.set(
            0,
            0,
            TracebackCell {
                score: 0,
                op: AlignmentOperation::Match(None),
            },
        );

        // construct the score matrix (O(n^2) space)
        let mut topo = Topo::new(&self.graph);
        while let Some(node) = topo.next(&self.graph) {
            // reference base and index
            let r = self.graph.raw_nodes()[node.index()].weight; // reference base at previous index
            let i = node.index() + 1;
            traceback.last = node;
            // iterate over the predecessors of this node
            let prevs: Vec<NodeIndex<usize>> =
                self.graph.neighbors_directed(node, Incoming).collect();
            // query base and its index in the DAG (traceback matrix rows)
            for (j_p, q) in query.iter().enumerate() {
                let j = j_p + 1;
                // match and deletion scores for the first reference base
                let max_cell = if prevs.is_empty() {
                    TracebackCell {
                        score: traceback.get(0, j - 1).score + self.scoring.match_fn.score(r, *q),
                        op: AlignmentOperation::Match(None),
                    }
                } else {
                    let mut max_cell = TracebackCell {
                        score: MIN_SCORE,
                        op: AlignmentOperation::Match(None),
                    };
                    for prev_node in &prevs {
                        let i_p: usize = prev_node.index() + 1; // index of previous node
                        let gap_penalty = self.determine_gap_penalty(traceback.get(i_p, j), AlignmentOperation::Del(Some((i_p - 1, i))));
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + self.scoring.match_fn.score(r, *q),
                                        op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + gap_penalty,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };

                let gap_penalty = self.determine_gap_penalty(traceback.get(i, j-1),AlignmentOperation::Ins(Some(i - 1)));
                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + gap_penalty,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
            }
        }

        traceback
    }
       


    /*/// Experimental: return sequence of traversed edges
    ///
    /// Only supports alignments for sequences that have already been added,
    /// so all operations must be Match.
    pub fn edges(&self, aln: Alignment) -> Vec<usize> {
        let mut path: Vec<usize> = vec![];
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut _i: usize = 0;
        for op in aln.operations {
            match op {
                AlignmentOperation::Match(None) => {
                    _i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(p);
                    let edge = self.graph.find_edge(prev, node).unwrap();
                    path.push(edge.index());
                    prev = NodeIndex::new(p);
                    _i += 1;
                }
                AlignmentOperation::Ins(None) => {}
                AlignmentOperation::Ins(Some(_)) => {}
                AlignmentOperation::Del(_) => {}
            }
        }
        path
    }*/

    /// return the consensus of the POA
    /// and also
    /// return the sequences associated to edges in the consensus
    pub fn consensus(&mut self, params : &Params) -> (Vec<u64>, Vec<String>) {
        let debug = params.debug; // prints POA graph in dot format

        let mut cns = Vec::new();
        let mut cns_edges_seqs = Vec::new();
        let nodes = self.consensus_path(&params).to_vec();
        for i in 0..(nodes.len())
        {
            cns.push(*self.graph.node_weight(nodes[i]).unwrap() as u64);
            if i < nodes.len()-1
            {
                let edge = self.graph.find_edge(nodes[i],nodes[i+1]).unwrap();
                cns_edges_seqs.push(self.graph.edge_weight(edge).unwrap().1.clone());
            }
        }
        //println!("[in poa] consensus()/consensus_edge_seqs() lens: {} / {}", cns.len(), cns_edges_seqs.len());

        if debug 
        {
            // debug: print graph in dot format
            let printable_graph = self.graph.filter_map(|_, w| Some(format!("{:?}", w)), |_, w| if w.0 > 1 { Some(format!("{:?}", w.0))} else {None});
            println!("{}", Dot::with_config(&printable_graph, &[]));
        }

        (cns.to_vec(), cns_edges_seqs.to_vec())
    }

    pub fn consensus_path(&mut self, params : &Params) -> Vec<NodeIndex<usize>> {
        // If we've added no new nodes or edges since the last call, sort first
        let node_idx = toposort(&self.graph, None).unwrap();
        // For each node find the best predecessor by edge-weight, breaking ties with path-weight
        let mut scores = HashMap::<NodeIndex<usize>, i32>::new();
        let mut next_in_path: HashMap<NodeIndex<usize>, Option<NodeIndex<usize>>> = HashMap::new();
        for node in node_idx.iter().rev() {
            let mut best_neighbor = None::<NodeIndex<usize>>;
            let mut best_weights = (0, 0); // (Edge-weight, Path-weight)
            for e_ref in self.graph.edges(*node) {
                let mut weight = e_ref.weight().0;
                if weight < params.t as i32 {
                    weight = 0;
                }
                let target = e_ref.target();

                if (weight, *scores.entry(target).or_insert(0)) > best_weights {
                    best_neighbor = Some(target);
                    best_weights = (weight, *scores.entry(target).or_insert(0));
                }
            }

            scores.insert(*node, best_weights.0 + best_weights.1);
            next_in_path.insert(*node, best_neighbor);
        }

        // Find the start position by locating the highest scoring node
        let mut start_node = None::<NodeIndex<usize>>;
        let mut best_score = 0;
        for (node, score) in &scores {
            if score > &best_score {
                start_node = Some(*node);
                best_score = *score;
            }
        }

        // Traverse the graph from the start node, recording the path of nodes taken
        let mut cns_path = Vec::new();
        let mut curr_node = start_node;
        loop {
            if curr_node.is_some() {
                let node = curr_node.unwrap();
                cns_path.push(node);
                curr_node = *next_in_path.get(&node).unwrap();
            } else {
                break;
            }
        }

        cns_path
    }
    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    ///
    pub fn add_alignment(&mut self, aln: &Alignment, seq: Vec<u64>, seq_str: Option<&String>, seq_minim_pos: Option<&Vec<usize>>) {
        let debug = false;
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut prev_i = 0;
        let mut i: usize = aln.ystart;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    let between_minims = String::from(&seq_str.unwrap()[seq_minim_pos.unwrap()[prev_i] as usize..seq_minim_pos.unwrap()[i] as usize]);
                    if debug {
                        println!("Weight : {}", self.graph.raw_nodes()[*p].weight);
                        println!("Node : {:?}", seq[i]);
                    }
                    if seq[i] != self.graph.raw_nodes()[*p].weight {
                        if debug {println!("node added");}
                        let node = self.graph.add_node(seq[i]);
                        self.graph.add_edge(prev, node, (1, between_minims));
                        //let edge = self.graph.find_edge(prev, node).unwrap(); 
                        // Rayan: I don't see the point of that above line
                        prev = node;
                        prev_i = i;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                if debug {println!("edge incremented");}
                                self.graph.edge_weight_mut(edge).unwrap().0 += 1;
                            }
                            None => {
                                if debug {println!("edge added");}
                                // where the previous node was newly added
                                self.graph.add_edge(prev, node, (1, between_minims));
                                //let edge = self.graph.find_edge(prev, node).unwrap();
                                //println!("Edge added {:?}", (self.node_index[&prev], self.node_index[&node]));
                            }
                        }
                        prev = node;
                        prev_i = i;
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    i += 1;
                }
               AlignmentOperation::Ins(Some(_)) => {
                    let node = self.graph.add_node(seq[i]);
                    let between_minims = String::from(&seq_str.unwrap()[seq_minim_pos.unwrap()[prev_i] as usize..seq_minim_pos.unwrap()[i] as usize]);
                    self.graph.add_edge(prev, node, (1, between_minims));
                    //let edge = self.graph.find_edge(prev, node).unwrap();
                    prev = node;
                    prev_i = i;
                    i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }

    /// Return a pretty formatted alignment, but it is a very rough approximation (we return a linear
    /// representation but keep in mind the template is a graph)
    /// interpretation:
    /// M = match, X mismatch
    /// D = deletion, I = insertion, i = some other strange kind of insertion
    pub fn pretty(&self, aln: &Alignment, seq: &[u64]) -> String{
        let mut x_pretty = String::new();
        let mut y_pretty = String::new();
        let mut inb_pretty = String::new();

        if !aln.operations.is_empty() {
            let mut x_add = String::new();
            let mut y_add = String::new();
            let mut inb_add = String::new();

            let mut i: usize = aln.ystart;
            for op in aln.operations.iter() {
                match op {
                    AlignmentOperation::Match(None) => {
                        // I think it's the start of the alignment, as per the alignment() function
                        // that breaks here
                        i += 1;
                    },
                    AlignmentOperation::Match(Some((_, p))) => {
                        if seq[i] != self.graph.raw_nodes()[*p].weight {
                            x_add   = String::from("X");//x[i];
                            inb_add = String::from("_");
                            y_add   = String::from("X");//seq[i].to_string().clone();
                        }
                        else {
                            x_add   = String::from("M");//x[i];
                            inb_add = String::from("_");
                            y_add   = String::from("M");//seq[i].to_string().clone();
                        }
                        i += 1;
                    },
                    AlignmentOperation::Del(_)=> {
                        x_add   = String::from("-");
                        inb_add = String::from("_");
                        y_add   = String::from("D");//seq[i].to_string().clone();
                    },
                    AlignmentOperation::Ins(_) => {
                        x_add   = String::from("I");//x[i];
                        inb_add = String::from("_");
                        y_add   = String::from("-");
                        i += 1;
                    },
                    /*AlignmentOperation::Ins(Some(_)) => {
                        x_add   = String::from("i");//x[i];
                        inb_add = String::from("_");
                        y_add   = String::from("-");
                        i += 1;
                    }*/

                }

                x_pretty.push_str(&x_add);
                inb_pretty.push_str(&inb_add);
                y_pretty.push_str(&y_add);
            }
        }
        format!("{}\n{}\n{}", x_pretty, inb_pretty, y_pretty)

    }
}
