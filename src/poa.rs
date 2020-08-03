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
//! # Example
//!
//! ```
//! use bio::alignment::poa::*;
//! use bio::alignment::pairwise::Scoring;
//!
//! let x = b"AAAAAAA";
//! let y = b"AABBBAA";
//! let z = b"AABCBAA";
//!
//! let scoring = Scoring::new(-1, 0, |a: u8, b: u8| if a == b { 1i32 } else { -1i32 });
//! let mut aligner = Aligner::new(scoring, x);
//! // z differs from x in 3 locations
//! assert_eq!(aligner.global(z).alignment().score, 1);
//! aligner.global(y).add_to_graph();
//! // z differs from x and y's partial order alignment by 1 base
//! assert_eq!(aligner.global(z).alignment().score, 5);
//! ```
//!

use std::cmp::{max, Ordering};


use std::collections::HashMap;

use petgraph::graph::NodeIndex;
use petgraph::visit::Topo;

use petgraph::{Directed, Graph, Incoming, Outgoing};

pub type POAGraph = Graph<u32, i32, Directed, usize>;

// Unlike with a total order we may have arbitrary successors in the
// traceback matrix. I have not yet figured out what the best level of
// detail to store is, so Match and Del operations remember In and Out
// nodes on the reference graph.
pub const MIN_SCORE: i32 = -858_993_459;

/// Trait required to instantiate a Scoring instance
pub trait MatchFunc {
    fn score(&self, a: u32, b: u32) -> i32;
}

/// A concrete data structure which implements trait MatchFunc with constant
/// match and mismatch scores
#[derive(Debug, Clone)]
pub struct MatchParams {
    pub match_score: i32,
    pub mismatch_score: i32,
}

impl MatchParams {
    /// Create new MatchParams instance with given match and mismatch scores
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
    }
}

impl MatchFunc for MatchParams {
    #[inline]
    fn score(&self, a: u32, b: u32) -> i32 {
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
    F: Fn(u32, u32) -> i32,
{
    fn score(&self, a: u32, b: u32) -> i32 {
        (self)(a, b)
    }
}

/// Details of scoring are encapsulated in this structure.
/// An affine gap score model is used so that the gap score for a length 'k' is:
/// GapScore(k) = gap_open + gap_extend * k
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
    /// Create new Scoring instance with given gap open, gap extend penalties
    /// match and mismatch scores. The clip penalties are set to MIN_SCORE by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_score` - the score for a match
    /// * `mismatch_score` - the score for a mismatch
    pub fn from_scores(
        gap_open: i32,
        gap_extend: i32,
        match_score: i32,
        mismatch_score: i32,
    ) -> Self {
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
    }
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

    /// Sets the prefix and suffix clipping penalties for x to the input value
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
    }
}

#[derive(Debug, Clone)]
pub enum AlignmentOperation {
    Match(Option<(usize, usize)>),
    Del(Option<(usize, usize)>),
    Ins(Option<usize>),
}

pub struct Alignment {
    pub score: i32,
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

    /// Populate the edges of the traceback matrix
    fn initialize_scores(&mut self, gap_open: i32) {
        for (i, row) in self
            .matrix
            .iter_mut()
            .enumerate()
            .take(self.rows + 1)
            .skip(1)
        {
            // TODO: these should be -1 * distance from head node
            row[0] = TracebackCell {
                score: (i as i32) * gap_open, // gap_open penalty
                op: AlignmentOperation::Del(None),
            };
        }
        for j in 1..=self.cols {
            self.matrix[0][j] = TracebackCell {
                score: (j as i32) * gap_open,
                op: AlignmentOperation::Ins(None),
            };
        }
    }

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

    pub fn print(&self, g: &Graph<u32, i32, Directed, usize>, query: Vec<u32>) {
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
    }
    

    pub fn alignment(&self) -> Alignment {
        // optimal AlignmentOperation path
        let mut ops: Vec<AlignmentOperation> = vec![];

        // Now backtrack through the matrix to construct an optimal path
        let mut i = self.last.index() + 1;
        let mut j = self.cols;

        while i > 0 && j > 0 {
            // push operation and edge corresponding to (one of the) optimal
            // routes
            ops.push(self.matrix[i][j].op.clone());
            match self.matrix[i][j].op {
                AlignmentOperation::Match(Some((p, _))) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Del(Some((p, _))) => {
                    i = p + 1;
                }
                AlignmentOperation::Ins(Some(p)) => {
                    i = p + 1;
                    j -= 1;
                }
                AlignmentOperation::Match(None) => {
                    break;
                }
                AlignmentOperation::Del(None) => {
                    j -= 1;
                }
                AlignmentOperation::Ins(None) => {
                    i -= 1;
                }
            }
        }

        ops.reverse();

        Alignment {
            score: self.matrix[self.last.index() + 1][self.cols].score,
            operations: ops,
        }
    }
}

/// A partially ordered aligner builder
///
/// Uses consuming builder pattern for constructing partial order alignments with method chaining
pub struct Aligner<F: MatchFunc> {
    pub traceback: Traceback,
    query: Vec<u32>,
    pub poa: Poa<F>,
}

impl<F: MatchFunc> Aligner<F> {
    /// Create new instance.
    pub fn new(scoring: Scoring<F>, reference: Vec<u32>) -> Self {
        Aligner {
            traceback: Traceback::new(),
            query: reference.to_vec(),
            poa: Poa::from_string(scoring, reference),
        }
    }

    /// Add the alignment of the last query to the graph.
    pub fn add_to_graph(&mut self) -> &mut Self {
        let alignment = self.traceback.alignment();
        self.poa.add_alignment(&alignment, self.query.to_vec());
        self
    }

    /// Return alignment of last added query against the graph.
    pub fn alignment(&self) -> Alignment {
        self.traceback.alignment()
    }

    /// Globally align a given query against the graph.
    pub fn global(&mut self, query: Vec<u32>) -> &mut Self {
        self.query = query.to_vec();
        self.traceback = self.poa.global(query);
        self
    }

    /// Return alignment graph.
    pub fn graph(&self) -> &POAGraph {
        &self.poa.graph
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
    pub node_index : HashMap<NodeIndex<usize>, u32>,
    pub edge_index : HashMap<usize, (u32, u32)>,
}

impl<F: MatchFunc> Poa<F> {
    /// Create a new aligner instance from the directed acyclic graph of another.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `poa` - the partially ordered reference alignment
    ///
    pub fn new(scoring: Scoring<F>, graph: POAGraph) -> Self {
        let mut node_index : HashMap<NodeIndex<usize>, u32> = HashMap::new();
        let mut edge_index : HashMap<usize, (u32, u32)> = HashMap::new();

        Poa { scoring, graph, node_index, edge_index }
    }

    /// Create a new POA graph from an initial reference sequence and alignment penalties.
    ///
    /// # Arguments
    ///
    /// * `scoring` - the score struct
    /// * `reference` - a reference TextSlice to populate the initial reference graph
    ///
    pub fn from_string(scoring: Scoring<F>, seq: Vec<u32>) -> Self {
        let mut node_index : HashMap<NodeIndex<usize>, u32> = HashMap::new();
        let mut edge_index : HashMap<usize, (u32, u32)> = HashMap::new();
        let mut graph: Graph<u32, i32, Directed, usize> =
            Graph::with_capacity(seq.len(), seq.len() - 1);
        let mut prev: NodeIndex<usize> = graph.add_node(seq[0]);
        node_index.insert(prev, seq[0]);
        let mut node: NodeIndex<usize>;
        for base in seq.iter().skip(1) {
            node = graph.add_node(*base);
            node_index.insert(node, *base);
            graph.add_edge(prev, node, 1);
            let edge = graph.find_edge(prev, node).unwrap();
            edge_index.insert(edge.index(),(node_index[&prev], node_index[&node]));
            prev = node;
        }

        Poa { scoring, graph, node_index, edge_index }
    }

    /// A global Needleman-Wunsch aligner on partially ordered graphs.
    ///
    /// # Arguments
    /// * `query` - the query TextSlice to align against the internal graph member
    ///

    pub fn global(&self, query: Vec<u32>) -> Traceback {
        assert!(self.graph.node_count() != 0);

        // dimensions of the traceback matrix
        let (m, n) = (self.graph.node_count(), query.len());
        let mut traceback = Traceback::with_capacity(m, n);
        traceback.initialize_scores(self.scoring.gap_open);

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
                        max_cell = max(
                            max_cell,
                            max(
                                TracebackCell {
                                    score: traceback.get(i_p, j - 1).score
                                        + self.scoring.match_fn.score(r, *q),
                                    op: AlignmentOperation::Match(Some((i_p - 1, i - 1))),
                                },
                                TracebackCell {
                                    score: traceback.get(i_p, j).score + self.scoring.gap_open,
                                    op: AlignmentOperation::Del(Some((i_p - 1, i))),
                                },
                            ),
                        );
                    }
                    max_cell
                };

                let score = max(
                    max_cell,
                    TracebackCell {
                        score: traceback.get(i, j - 1).score + self.scoring.gap_open,
                        op: AlignmentOperation::Ins(Some(i - 1)),
                    },
                );
                traceback.set(i, j, score);
            }
        }

        traceback
    }

    /// Experimental: return sequence of traversed edges
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
    }

    /// Incorporate a new sequence into a graph from an alignment
    ///
    /// # Arguments
    ///
    /// * `aln` - The alignment of the new sequence to the graph
    /// * `seq` - The sequence being incorporated
    ///
    pub fn find_consensus(&mut self, cons_scores : Vec<i32>, cons_next : Vec<i32>) -> Vec<u32> {
        let mut consensus = Vec::<u32>::new();
        let max = cons_scores.iter().max().unwrap();
        let mut pos = cons_scores.iter().position(|element| element == max).unwrap() as usize;
        let mut first : bool = true;
        while (pos != 0 && !first) || first {
            let min = self.graph.raw_nodes()[pos].weight;
            //println!("Index : {}, min : {}", pos, min);
            consensus.push(min);
            first = false;
            pos = cons_next[pos] as usize;
        }
        consensus

    }
    pub fn set_max_weight_scores(&mut self) -> (Vec<i32>, Vec<i32>) {
        let mut cons_scores = Vec::<i32>::new();
        let mut cons_next = Vec::<i32>::new();
        for i in 0..self.graph.node_count() {
            cons_scores.push(-1);
            cons_next.push(-1);
        }
        for i in 0..self.graph.node_count() {
            let prev_index = self.graph.node_count()-i-1;
            let prev = NodeIndex::new(self.graph.node_count()-i-1);
            let next: Vec<NodeIndex<usize>> = self.graph.neighbors_directed(prev, Outgoing).collect();
            let mut max_weight = -1;
            let mut max_node = NodeIndex::new(0 as usize);
            let mut count = 0;
            let mut max_weight_nodes = Vec::<NodeIndex<usize>>::new();
            for node in next.iter() {
                let edge = self.graph.find_edge(prev, *node).unwrap();
                let weight = self.graph.edge_weight_mut(edge).unwrap();
                if *weight == max_weight {
                    max_node = node.clone();
                    max_weight_nodes.push(max_node);
                }
                else if *weight > max_weight {
                    max_weight = *weight;
                    max_weight_nodes = Vec::<NodeIndex<usize>>::new();
                    max_node = node.clone();
                    max_weight_nodes.push(max_node);
                }
            }
            if max_weight_nodes.len() != 1 {
                //println!("Here");
                let mut max_score = 0;
                let mut max_score_node = NodeIndex::new(0 as usize);
                for max_node in max_weight_nodes.iter() {
                    let max_index = max_node.index();
                    let score = &cons_scores[max_index];
                    if score > &max_score {
                        max_score = *score;
                        max_score_node = max_node.clone();
                    }
                }
                let final_index = max_score_node.index();
                let final_score = max_weight + max_score;
                cons_scores[prev_index] = final_score;
                cons_next[prev_index] = final_index as i32;
            }
            else {
                let final_index = max_node.index();
                let final_score = max_weight + cons_scores[final_index];
                cons_scores[prev_index] = final_score;
                cons_next[prev_index] = final_index as i32;
            }
            
            //println!("Index : {}, score : {}, pointing to : {}", prev_index, cons_scores[prev_index], cons_next[prev_index]);
        }
        (cons_scores, cons_next)


    }
    pub fn add_alignment(&mut self, aln: &Alignment, seq: Vec<u32>) {
        let mut prev: NodeIndex<usize> = NodeIndex::new(0);
        let mut i: usize = 0;
        for op in aln.operations.iter() {
            match op {
                AlignmentOperation::Match(None) => {
                    i += 1;
                }
                AlignmentOperation::Match(Some((_, p))) => {
                    let node = NodeIndex::new(*p);
                    //println!("Weight : {}", self.graph.raw_nodes()[*p].weight);
                    //println!("Node : {:?}", seq[i]);
                    if (seq[i] != self.graph.raw_nodes()[*p].weight) {
                        let node = self.graph.add_node(seq[i]);
                        self.node_index.insert(node, seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        let edge = self.graph.find_edge(prev, node).unwrap();
                        self.edge_index.insert(edge.index(),(self.node_index[&prev], self.node_index[&node]));
                        prev = node;
                    } else {
                        // increment node weight
                        match self.graph.find_edge(prev, node) {
                            Some(edge) => {
                                *self.graph.edge_weight_mut(edge).unwrap() += 1;
                            }
                            None => {
                                // where the previous node was newly added
                                self.graph.add_edge(prev, node, 1);
                                let edge = self.graph.find_edge(prev, node).unwrap();
                                //println!("Edge added {:?}", (self.node_index[&prev], self.node_index[&node]));
                                self.edge_index.insert(edge.index(),(self.node_index[&prev], self.node_index[&node]));
                            }
                        }
                        let node = NodeIndex::new(*p);
                        self.node_index.insert(node, seq[i]);
                        prev = node;
                    }
                    i += 1;
                }
                AlignmentOperation::Ins(None) => {
                    i += 1;
                }
                AlignmentOperation::Ins(Some(_)) => {
                        let node = self.graph.add_node(seq[i]);
                        self.node_index.insert(node, seq[i]);
                        self.graph.add_edge(prev, node, 1);
                        let edge = self.graph.find_edge(prev, node).unwrap();
                        self.edge_index.insert(edge.index(),(self.node_index[&prev], self.node_index[&node]));
                        prev = node;
                        i += 1;
                }
                AlignmentOperation::Del(_) => {} // we should only have to skip over deleted nodes
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::pairwise::Scoring;
    use petgraph::graph::NodeIndex;

    #[test]
    fn test_init_graph() {
        // sanity check for String -> Graph

        let scoring = Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });
        let poa = Poa::from_string(scoring, b"123456789");
        assert!(poa.graph.is_directed());
        assert_eq!(poa.graph.node_count(), 9);
        assert_eq!(poa.graph.edge_count(), 8);
    }

    #[test]
    fn test_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });
        // examples from the POA paper
        //let _seq1 = b"PKMIVRPQKNETV";
        //let _seq2 = b"THKMLVRNETIM";
        let poa = Poa::from_string(scoring, b"GATTACA");
        let alignment = poa.global(b"GCATGCU").alignment();
        assert_eq!(alignment.score, 0);

        let alignment = poa.global(b"GCATGCUx").alignment();
        assert_eq!(alignment.score, -1);

        let alignment = poa.global(b"xCATGCU").alignment();
        assert_eq!(alignment.score, -2);
    }

    #[test]
    fn test_branched_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });
        let seq1 = b"TTTTT";
        let seq2 = b"TTATT";
        let mut poa = Poa::from_string(scoring, seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.global(seq2).alignment();
        assert_eq!(alignment.score, 3);
    }

    #[test]
    fn test_alt_branched_alignment() {
        let scoring = Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });

        let seq1 = b"TTCCTTAA";
        let seq2 = b"TTTTGGAA";
        let mut poa = Poa::from_string(scoring, seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'A');
        let node2 = poa.graph.add_node(b'A');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, tail, 1);
        let alignment = poa.global(seq2).alignment();
        poa.add_alignment(&alignment, seq2);
        assert_eq!(poa.graph.edge_count(), 14);
        assert!(poa
            .graph
            .contains_edge(NodeIndex::new(5), NodeIndex::new(10)));
        assert!(poa
            .graph
            .contains_edge(NodeIndex::new(11), NodeIndex::new(6)));
    }

    #[test]
    fn test_insertion_on_branch() {
        let scoring = Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });

        let seq1 = b"TTCCGGTTTAA";
        let seq2 = b"TTGGTATGGGAA";
        let seq3 = b"TTGGTTTGCGAA";
        let mut poa = Poa::from_string(scoring, seq1);
        let head: NodeIndex<usize> = NodeIndex::new(1);
        let tail: NodeIndex<usize> = NodeIndex::new(2);
        let node1 = poa.graph.add_node(b'C');
        let node2 = poa.graph.add_node(b'C');
        let node3 = poa.graph.add_node(b'C');
        poa.graph.add_edge(head, node1, 1);
        poa.graph.add_edge(node1, node2, 1);
        poa.graph.add_edge(node2, node3, 1);
        poa.graph.add_edge(node3, tail, 1);
        let alignment = poa.global(seq2).alignment();
        assert_eq!(alignment.score, 2);
        poa.add_alignment(&alignment, seq2);
        let alignment2 = poa.global(seq3).alignment();

        assert_eq!(alignment2.score, 10);
    }

    #[test]
    fn test_poa_method_chaining() {
        let scoring = Scoring::new(-1, 0, |a: u32, b: u32| if a == b { 1i32 } else { -1i32 });
        let mut aligner = Aligner::new(scoring, b"TTCCGGTTTAA");
        aligner
            .global(b"TTGGTATGGGAA")
            .add_to_graph()
            .global(b"TTGGTTTGCGAA")
            .add_to_graph();
        assert_eq!(aligner.alignment().score, 10);
    }
}