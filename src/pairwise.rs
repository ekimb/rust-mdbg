// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use std::cmp::max;
use std::i32;
use std::iter::repeat;

// Types for representing pairwise sequence alignments

pub type TextSlice<'a> = &'a [u64];

/// Alignment operations supported are match, substitution, insertion, deletion
/// and clipping. Clipping is a special boundary condition where you are allowed
/// to clip off the beginning/end of the sequence for a fixed clip penalty. The
/// clip penalty could be different for the two sequences x and y, and the
/// clipping operations on both are distinguishable (Xclip and Yclip). The usize
/// value associated with the clipping operations are the lengths clipped. In case
/// of standard modes like Global, Semi-Global and Local alignment, the clip operations
/// are filtered out

#[derive(Eq, PartialEq, Debug, Copy, Clone)]//, Serialize, Deserialize)]
pub enum AlignmentOperation {
    Match,
    Subst,
    Del,
    Ins,
    Xclip(usize),
    Yclip(usize),
}

/// The modes of alignment supported by the aligner include standard modes such as
/// Global, Semi-Global and Local alignment. In addition to this, user can also invoke
/// the custom mode. In the custom mode, users can explicitly specify the clipping penalties
/// for prefix and suffix of strings 'x' and 'y' independently. Under the hood the standard
/// modes are implemented as special cases of the custom mode with the clipping penalties
/// appropriately set.
///
/// The default alignment mode is Global.

#[derive(Debug, PartialEq, Eq, Copy, Clone)]//, Serialize, Deserialize)]
pub enum AlignmentMode {
    //Local,
    Semiglobal,
    Global,
    Custom,
}

impl Default for AlignmentMode {
    fn default() -> Self {
        AlignmentMode::Global
    }
}

/// We consider alignment between two sequences x and  y. x is the query or read sequence
/// and y is the reference or template sequence. An alignment, consisting of a score,
/// the start and end position of the alignment on sequence x and sequence y, the
/// lengths of sequences x and y, and the alignment edit operations. The start position
/// and end position of the alignment does not include the clipped regions. The length
/// of clipped regions are already encapsulated in the Alignment Operation.

#[derive(Debug, Eq, PartialEq, Clone)]//Serialize, Deserialize, Default)]
pub struct Alignment {
    /// Smith-Waterman alignment score
    pub score: i32,
    /// Start position of alignment in reference
    pub ystart: usize,
    /// Start position of alignment in query
    pub xstart: usize,
    /// End position of alignment in reference
    pub yend: usize,
    /// End position of alignment in query
    pub xend: usize,
    /// Length of the reference sequence
    pub ylen: usize,
    /// Length of the query sequence
    pub xlen: usize,
    /// Vector of alignment operations
    pub operations: Vec<AlignmentOperation>,
    pub mode: AlignmentMode,
}

impl Alignment {
    /// Calculate the cigar string from the alignment struct. x is the target string
    ///
    /// # Example
    ///
    /// ```
    /// use bio_types::alignment::{Alignment,AlignmentMode};
    /// use bio_types::alignment::AlignmentOperation::{Match, Subst, Ins, Del};
    /// let alignment = Alignment {
    ///     score: 5,
    ///     xstart: 3,
    ///     ystart: 0,
    ///     xend: 9,
    ///     yend: 10,
    ///     ylen: 10,
    ///     xlen: 10,
    ///     operations: vec![Match, Match, Match, Subst, Ins, Ins, Del, Del],
    ///     mode: AlignmentMode::Semiglobal
    /// };
    /// assert_eq!(alignment.cigar(false), "3S3=1X2I2D1S");
    /// ```
    /*pub fn cigar(&self, hard_clip: bool) -> String {
        match self.mode {
            AlignmentMode::Global => panic!("CIGAR function not supported for Global Alignment mode."),
            AlignmentMode::Local => panic!("CIGAR function not supported for Local Alignment mode."),
            _ => {}
        }
        let clip_str = if hard_clip {"H"} else {"S"};
        let add_op = |op: AlignmentOperation, k, cigar: &mut String| match op {
            AlignmentOperation::Match => cigar.push_str(&format!("{}{}", k, "=")),
            AlignmentOperation::Subst => cigar.push_str(&format!("{}{}", k, "X")),
            AlignmentOperation::Del => cigar.push_str(&format!("{}{}", k, "D")),
            AlignmentOperation::Ins => cigar.push_str(&format!("{}{}", k, "I")),
            _ => {}
        };
        let mut cigar = "".to_owned();
        if self.operations.is_empty() {
            return cigar;
        }
        let mut last = self.operations[0];
        if self.xstart > 0 {
            cigar.push_str(&format!("{}{}", self.xstart, clip_str))
        }
        let mut k = 1;
        for &op in self.operations[1..].iter() {
            if op == last {
                k += 1;
            } else {
                add_op(last, k, &mut cigar);
                k = 1;
            }
            last = op;
        }
        add_op(last, k, &mut cigar);
        if self.xlen > self.xend {
            cigar.push_str(&format!("{}{}", self.xlen - self.xend, clip_str))
        }
        cigar
    }

    /// Returns the optimal path in the alignment matrix
    pub fn path(&self) -> Vec<(usize, usize, AlignmentOperation)> {
        let mut path = Vec::new();

        if !self.operations.is_empty() {
            let last = match self.mode {
                AlignmentMode::Custom => (self.xlen, self.ylen),
                _ => (self.xend, self.yend),
            };
            let mut x_i = last.0;
            let mut y_i = last.1;

            let mut ops = self.operations.clone();
            ops.reverse();

            // Process the alignment.
            for i in ops {
                path.push((x_i, y_i, i));
                match i {
                    AlignmentOperation::Match => {
                        x_i -= 1;
                        y_i -= 1;
                    }
                    AlignmentOperation::Subst => {
                        x_i -= 1;
                        y_i -= 1;
                    }
                    AlignmentOperation::Del => {
                        y_i -= 1;
                    }
                    AlignmentOperation::Ins => {
                        x_i -= 1;
                    }
                    AlignmentOperation::Xclip(len) => {
                        x_i -= len;
                    }
                    AlignmentOperation::Yclip(len) => {
                        y_i -= len;
                    }
                }
            }
        }
        path.reverse();
        path
    }*/

    /// Filter out Xclip and Yclip operations from the list of operations. Useful
    /// when invoking the standard modes.
    pub fn filter_clip_operations(&mut self) {
        use self::AlignmentOperation::{Del, Ins, Match, Subst};
        self.operations
            .retain(|x| (*x == Match || *x == Subst || *x == Ins || *x == Del));
    }

    /*/// Number of bases in reference sequence that are aligned
    pub fn y_aln_len(&self) -> usize {
        self.yend - self.ystart
    }

    /// Number of bases in query sequence that are aigned
    pub fn x_aln_len(&self) -> usize {
        self.xend - self.xstart
    }*/
}


//// --- snip for bio-types


/// Value to use as a 'negative infinity' score. Should be close to `i32::MIN`,
/// but avoid underflow when used with reasonable scoring parameters or even
/// adding two negative infinities. Use ~ `0.4 * i32::MIN`
pub const MIN_SCORE: i32 = -858_993_459;

/// Trait required to instantiate a Scoring instance
pub trait MatchFunc {
    fn score(&self, a: u64, b: u64) -> i32;
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
    fn score(&self, a: u64, b: u64) -> i32 {
        if a == b {
            self.match_score
        } else {
            self.mismatch_score
        }
    }
}

/// The trait Matchfunc is also implemented for Fn(u64, u64) -> i32 so that Scoring
/// can be instantiated using closures and custom user defined functions
impl<F> MatchFunc for F
where
    F: Fn(u64, u64) -> i32,
{
    fn score(&self, a: u64, b: u64) -> i32 {
        (self)(a, b)
    }
}

/// Details of scoring are encapsulated in this structure.
///
/// An [affine gap score model](https://en.wikipedia.org/wiki/Gap_penalty#Affine)
/// is used so that the gap score for a length `k` is:
/// `GapScore(k) = gap_open + gap_extend * k`
#[derive(Debug, Clone)]
pub struct Scoring<F: MatchFunc> {
    pub gap_open: i32,
    pub gap_extend: i32,
    pub match_fn: F,
#[allow(dead_code)]
    pub match_scores: Option<(i32, i32)>,
    pub xclip_prefix: i32,
    pub xclip_suffix: i32,
    pub yclip_prefix: i32,
    pub yclip_suffix: i32,
}

impl Scoring<MatchParams> {
    /*/// Create new Scoring instance with given gap open, gap extend penalties
    /// match and mismatch scores. The clip penalties are set to `MIN_SCORE` by default
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
    /// and the score function. The clip penalties are set to [`MIN_SCORE`](constant.MIN_SCORE.html) by default
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should not be positive)
    /// * `gap_extend` - the score for extending a gap (should not be positive)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
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
    }*/

    /*/// Sets the prefix clipping penalty for x to the input value
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
    }*/

    /*/// Sets the suffix clipping penalty for x to the input value
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
    }*/

    /*/// Sets the prefix and suffix clipping penalties for y to the input value
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
    /// ```*/
    /*pub fn yclip(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self.yclip_suffix = penalty;
        self
    }*/

    /*/// Sets the prefix clipping penalty for y to the input value
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
    /// ```*/
    /*pub fn yclip_prefix(mut self, penalty: i32) -> Self {
        assert!(penalty <= 0, "Clipping penalty can't be positive");
        self.yclip_prefix = penalty;
        self
    }*/

    /*/// Sets the suffix clipping penalty for y to the input value
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

/// A generalized Smith-Waterman aligner.
///
/// `M(i,j)` is the best score such that `x[i]` and `y[j]` ends in a match (or substitution)
/// ```ignore
///              .... A   G  x_i
///              .... C   G  y_j
/// ```
/// `I(i,j)` is the best score such that `x[i]` is aligned with a gap
/// ```ignore
///              .... A   G  x_i
///              .... G  y_j  -
/// ```
/// This is interpreted as an insertion into `x` w.r.t reference `y`
///
/// `D(i,j)` is the best score such that `y[j]` is aligned with a gap
/// ```ignore
///              .... A  x_i  -
///              .... G   G  y_j
/// ```
/// This is interpreted as a deletion from `x` w.r.t reference `y`
///
/// `S(i,j)` is the best score for prefixes `x[0..i]`, `y[0..j]`
///
/// To save space, only two columns of these matrices are stored at
/// any point - the current column and the previous one. Moreover
/// `M(i,j)` is not explicitly stored
///
/// `Lx` is the optimal x suffix clipping lengths from each position of the
/// sequence y
///
/// `Ly` is the optimal y suffix clipping lengths from each position of the
/// sequence x
///
/// `Sn` is the last column of the matrix. This is needed to keep track of
/// suffix clipping scores
///
/// `traceback` - see [`bio::alignment::pairwise::TracebackCell`](struct.TracebackCell.html)
///
/// `scoring` - see [`bio::alignment::pairwise::Scoring`](struct.Scoring.html)
#[allow(non_snake_case)]
pub struct Aligner<F: MatchFunc> {
    I: [Vec<i32>; 2],
    D: [Vec<i32>; 2],
    S: [Vec<i32>; 2],
    Lx: Vec<usize>,
    Ly: Vec<usize>,
    Sn: Vec<i32>,
    traceback: Traceback,
    scoring: Scoring<F>,
}

/*fn seq_to_vec(seq: &[u8]) -> Vec<u64> {
    let mut res = Vec::new();
    for c in seq.iter()
    {
        res.push(match *c as char { 'A' => 1, 'C' => 2, 'G' => 3, 'T' => 4, 'U' => 5, 'x' => 6,_ => 0});
    }
    res
}*/

//const DEFAULT_ALIGNER_CAPACITY: usize = 200;

impl<F: MatchFunc> Aligner<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    /*pub fn new(gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        Aligner::with_capacity(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            gap_open,
            gap_extend,
            match_fn,
        )
    }*/

    /// Create new aligner instance. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions
    ///    (see also [`bio::alignment::pairwise::Scoring`](struct.Scoring.html))
    pub fn with_capacity(m: usize, n: usize, gap_open: i32, gap_extend: i32, match_fn: F) -> Self {
        assert!(gap_open <= 0, "gap_open can't be positive");
        assert!(gap_extend <= 0, "gap_extend can't be positive");

        Aligner {
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring: Scoring::new(gap_open, gap_extend, match_fn),
        }
    }

    /// Create new aligner instance with given the scoring struct
    ///
    /// # Arguments
    ///
    /// * `scoring` - the scoring struct (see bio::alignment::pairwise::Scoring)
    /*pub fn with_scoring(scoring: Scoring<F>) -> Self {
        Aligner::with_capacity_and_scoring(
            DEFAULT_ALIGNER_CAPACITY,
            DEFAULT_ALIGNER_CAPACITY,
            scoring,
        )
    }*/

    /// Create new aligner instance with scoring and size hint. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `x_exp` - the expected size of x
    /// * `y_exp` - the expected size of y
    /// * `scoring` - the scoring struct
    /*pub fn with_capacity_and_scoring(m: usize, n: usize, scoring: Scoring<F>) -> Self {
        assert!(scoring.gap_open <= 0, "gap_open can't be positive");
        assert!(scoring.gap_extend <= 0, "gap_extend can't be positive");
        assert!(
            scoring.xclip_prefix <= 0,
            "Clipping penalty (x prefix) can't be positive"
        );
        assert!(
            scoring.xclip_suffix <= 0,
            "Clipping penalty (x suffix) can't be positive"
        );
        assert!(
            scoring.yclip_prefix <= 0,
            "Clipping penalty (y prefix) can't be positive"
        );
        assert!(
            scoring.yclip_suffix <= 0,
            "Clipping penalty (y suffix) can't be positive"
        );

        Aligner {
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring,
        }
    }*/

    /// The core function to compute the alignment
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    pub fn custom(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        let (x_exp, y_exp) = (x.len(), y.len());
        self.traceback.init(x_exp, y_exp);

        // Set the initial conditions
        // We are repeating some work, but that's okay!
        for k in 0..2 {
            self.I[k].clear();
            self.D[k].clear();
            self.S[k].clear();

            self.D[k].extend(repeat(MIN_SCORE).take(x_exp + 1));
            self.I[k].extend(repeat(MIN_SCORE).take(x_exp + 1));
            self.S[k].extend(repeat(MIN_SCORE).take(x_exp + 1));

            self.S[k][0] = 0;

            if k == 0 {
                let mut tb = TracebackCell::new();
                tb.set_all(TB_START);
                self.traceback.set(0, 0, tb);
                self.Lx.clear();
                self.Lx.extend(repeat(0usize).take(y_exp + 1));
                self.Ly.clear();
                self.Ly.extend(repeat(0usize).take(x_exp + 1));
                self.Sn.clear();
                self.Sn.extend(repeat(MIN_SCORE).take(x_exp + 1));
                self.Sn[0] = self.scoring.yclip_suffix;
                self.Ly[0] = y_exp;
            }

            for i in 1..=x_exp {
                let mut tb = TracebackCell::new();
                tb.set_all(TB_START);
                if i == 1 {
                    self.I[k][i] = self.scoring.gap_open + self.scoring.gap_extend;
                    tb.set_i_bits(TB_START);
                } else {
                    // Insert all i characters
                    let i_score = self.scoring.gap_open + self.scoring.gap_extend * (i as i32);
                    let c_score =
                        self.scoring.xclip_prefix + self.scoring.gap_open + self.scoring.gap_extend; // Clip then insert
                    if i_score > c_score {
                        self.I[k][i] = i_score;
                        tb.set_i_bits(TB_INS);
                    } else {
                        self.I[k][i] = c_score;
                        tb.set_i_bits(TB_XCLIP_PREFIX);
                    }
                }

                if i == x_exp {
                    tb.set_s_bits(TB_XCLIP_SUFFIX);
                } else {
                    self.S[k][i] = MIN_SCORE;
                }

                if self.I[k][i] > self.S[k][i] {
                    self.S[k][i] = self.I[k][i];
                    tb.set_s_bits(TB_INS);
                }

                if self.scoring.xclip_prefix > self.S[k][i] {
                    self.S[k][i] = self.scoring.xclip_prefix;
                    tb.set_s_bits(TB_XCLIP_PREFIX);
                }

                // Track the score if we do a suffix clip (x) after this character
                if i != x_exp && self.S[k][i] + self.scoring.xclip_suffix > self.S[k][x_exp] {
                    self.S[k][x_exp] = self.S[k][i] + self.scoring.xclip_suffix;
                    self.Lx[0] = x_exp - i;
                }

                if k == 0 {
                    self.traceback.set(i, 0, tb);
                }
                // Track the score if we do suffix clip (y) from here
                if self.S[k][i] + self.scoring.yclip_suffix > self.Sn[i] {
                    self.Sn[i] = self.S[k][i] + self.scoring.yclip_suffix;
                    self.Ly[i] = y_exp;
                }
            }
        }

        for j in 1..=y_exp {
            let curr = j % 2;
            let prev = 1 - curr;

            {
                // Handle i = 0 case
                let mut tb = TracebackCell::new();
                self.I[curr][0] = MIN_SCORE;

                if j == 1 {
                    self.D[curr][0] = self.scoring.gap_open + self.scoring.gap_extend;
                    tb.set_d_bits(TB_START);
                } else {
                    // Delete all j characters
                    let d_score = self.scoring.gap_open + self.scoring.gap_extend * (j as i32);
                    let c_score =
                        self.scoring.yclip_prefix + self.scoring.gap_open + self.scoring.gap_extend;
                    if d_score > c_score {
                        self.D[curr][0] = d_score;
                        tb.set_d_bits(TB_DEL);
                    } else {
                        self.D[curr][0] = c_score;
                        tb.set_d_bits(TB_YCLIP_PREFIX);
                    }
                }
                if self.D[curr][0] > self.scoring.yclip_prefix {
                    self.S[curr][0] = self.D[curr][0];
                    tb.set_s_bits(TB_DEL);
                } else {
                    self.S[curr][0] = self.scoring.yclip_prefix;
                    tb.set_s_bits(TB_YCLIP_PREFIX);
                }

                if j == y_exp && self.Sn[0] > self.S[curr][0] {
                    // Check if the suffix clip score is better
                    self.S[curr][0] = self.Sn[0];
                    tb.set_s_bits(TB_YCLIP_SUFFIX);
                // Track the score if we do suffix clip (y) from here
                } else if self.S[curr][0] + self.scoring.yclip_suffix > self.Sn[0] {
                    self.Sn[0] = self.S[curr][0] + self.scoring.yclip_suffix;
                    self.Ly[0] = y_exp - j;
                }

                self.traceback.set(0, j, tb);
            }

            for i in 1..=x_exp {
                self.S[curr][i] = MIN_SCORE;
            }

            let q = y[j - 1];
            let xclip_score = self.scoring.xclip_prefix
                + max(
                    self.scoring.yclip_prefix,
                    self.scoring.gap_open + self.scoring.gap_extend * (j as i32),
                );
            for i in 1..x_exp + 1 {
                let p = x[i - 1];
                let mut tb = TracebackCell::new();

                let m_score = self.S[prev][i - 1] + self.scoring.match_fn.score(p, q);

                let i_score = self.I[curr][i - 1] + self.scoring.gap_extend;
                let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
                let best_i_score;
                if i_score > s_score {
                    best_i_score = i_score;
                    tb.set_i_bits(TB_INS);
                } else {
                    best_i_score = s_score;
                    tb.set_i_bits(self.traceback.get(i - 1, j).get_s_bits());
                }

                let d_score = self.D[prev][i] + self.scoring.gap_extend;
                let s_score = self.S[prev][i] + self.scoring.gap_open + self.scoring.gap_extend;
                let best_d_score;
                if d_score > s_score {
                    best_d_score = d_score;
                    tb.set_d_bits(TB_DEL);
                } else {
                    best_d_score = s_score;
                    tb.set_d_bits(self.traceback.get(i, j - 1).get_s_bits());
                }

                tb.set_s_bits(TB_XCLIP_SUFFIX);
                let mut best_s_score = self.S[curr][i];

                if m_score > best_s_score {
                    best_s_score = m_score;
                    tb.set_s_bits(if p == q { TB_MATCH } else { TB_SUBST });
                }

                if best_i_score > best_s_score {
                    best_s_score = best_i_score;
                    tb.set_s_bits(TB_INS);
                }

                if best_d_score > best_s_score {
                    best_s_score = best_d_score;
                    tb.set_s_bits(TB_DEL);
                }

                if xclip_score > best_s_score {
                    best_s_score = xclip_score;
                    tb.set_s_bits(TB_XCLIP_PREFIX);
                }

                let yclip_score = self.scoring.yclip_prefix
                    + self.scoring.gap_open
                    + self.scoring.gap_extend * (i as i32);
                if yclip_score > best_s_score {
                    best_s_score = yclip_score;
                    tb.set_s_bits(TB_YCLIP_PREFIX);
                }

                self.S[curr][i] = best_s_score;
                self.I[curr][i] = best_i_score;
                self.D[curr][i] = best_d_score;

                // Track the score if we do suffix clip (x) from here
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][x_exp] {
                    self.S[curr][x_exp] = self.S[curr][i] + self.scoring.xclip_suffix;
                    self.Lx[j] = x_exp - i;
                }

                // Track the score if we do suffix clip (y) from here
                if self.S[curr][i] + self.scoring.yclip_suffix > self.Sn[i] {
                    self.Sn[i] = self.S[curr][i] + self.scoring.yclip_suffix;
                    self.Ly[i] = y_exp - j;
                }

                self.traceback.set(i, j, tb);
            }
        }

        // Handle suffix clipping in the j=n case
        for i in 0..=x_exp {
            let j = y_exp;
            let curr = j % 2;
            if self.Sn[i] > self.S[curr][i] {
                self.S[curr][i] = self.Sn[i];
                self.traceback.get_mut(i, j).set_s_bits(TB_YCLIP_SUFFIX);
            }
            if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][x_exp] {
                self.S[curr][x_exp] = self.S[curr][i] + self.scoring.xclip_suffix;
                self.Lx[j] = x_exp - i;
                self.traceback.get_mut(x_exp, j).set_s_bits(TB_XCLIP_SUFFIX);
            }
        }

        // Since there could be a change in the last column of S,
        // recompute the last column of I as this could also change
        for i in 1..=x_exp {
            let j = y_exp;
            let curr = j % 2;
            let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            if s_score > self.I[curr][i] {
                self.I[curr][i] = s_score;
                let s_bit = self.traceback.get(i - 1, j).get_s_bits();
                self.traceback.get_mut(i, j).set_i_bits(s_bit);
            }
            if s_score > self.S[curr][i] {
                self.S[curr][i] = s_score;
                self.traceback.get_mut(i, j).set_s_bits(TB_INS);
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][x_exp] {
                    self.S[curr][x_exp] = self.S[curr][i] + self.scoring.xclip_suffix;
                    self.Lx[j] = x_exp - i;
                    self.traceback.get_mut(x_exp, j).set_s_bits(TB_XCLIP_SUFFIX);
                }
            }
        }

        let mut i = x_exp;
        let mut j = y_exp;
        let mut operations = Vec::with_capacity(x.len());
        let mut xstart: usize = 0usize;
        let mut ystart: usize = 0usize;
        let mut xend = x_exp;
        let mut yend = y_exp;

        let mut last_layer = self.traceback.get(i, j).get_s_bits();

        loop {
            let next_layer: u16;
            match last_layer {
                TB_START => break,
                TB_INS => {
                    operations.push(AlignmentOperation::Ins);
                    next_layer = self.traceback.get(i, j).get_i_bits();
                    i -= 1;
                }
                TB_DEL => {
                    operations.push(AlignmentOperation::Del);
                    next_layer = self.traceback.get(i, j).get_d_bits();
                    j -= 1;
                }
                TB_MATCH => {
                    operations.push(AlignmentOperation::Match);
                    next_layer = self.traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                TB_SUBST => {
                    operations.push(AlignmentOperation::Subst);
                    next_layer = self.traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                TB_XCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Xclip(i));
                    xstart = i;
                    i = 0;
                    next_layer = self.traceback.get(0, j).get_s_bits();
                }
                TB_XCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Xclip(self.Lx[j]));
                    i -= self.Lx[j];
                    xend = i;
                    next_layer = self.traceback.get(i, j).get_s_bits();
                }
                TB_YCLIP_PREFIX => {
                    operations.push(AlignmentOperation::Yclip(j));
                    ystart = j;
                    j = 0;
                    next_layer = self.traceback.get(i, 0).get_s_bits();
                }
                TB_YCLIP_SUFFIX => {
                    operations.push(AlignmentOperation::Yclip(self.Ly[i]));
                    j -= self.Ly[i];
                    yend = j;
                    next_layer = self.traceback.get(i, j).get_s_bits();
                }
                _ => panic!("Dint expect this!"),
            }
            last_layer = next_layer;
        }

        operations.reverse();
        Alignment {
            score: self.S[y_exp % 2][x_exp],
            ystart,
            xstart,
            yend,
            xend,
            ylen: y_exp,
            xlen: x_exp,
            operations,
            mode: AlignmentMode::Custom,
        }
    }

    /*pub fn global_str(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        self.global(&seq_to_vec(x),&seq_to_vec(y))
    }*/
    /// Calculate global alignment of x against y.
    /*pub fn global(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        // Store the current clip penalties
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = MIN_SCORE;
        self.scoring.yclip_suffix = MIN_SCORE;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Global;

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }*/

    /*pub fn semiglobal_str(&mut self, x: &[u8], y: &[u8]) -> Alignment {
        self.semiglobal(&seq_to_vec(x),&seq_to_vec(y))
    }*/

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        // Store the current clip penalties
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Semiglobal;

        // Filter out Xclip and Yclip from alignment.operations
        alignment.filter_clip_operations();

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    // Calculate local alignment of x against y.
    /*pub fn local(&mut self, x: TextSlice<'_>, y: TextSlice<'_>) -> Alignment {
        // Store the current clip penalties
        let clip_penalties = [
            self.scoring.xclip_prefix,
            self.scoring.xclip_suffix,
            self.scoring.yclip_prefix,
            self.scoring.yclip_suffix,
        ];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = 0;
        self.scoring.xclip_suffix = 0;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Local;

        // Filter out Xclip and Yclip from alignment.operations
        alignment.filter_clip_operations();

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }*/
}

/// Packed representation of one cell of a Smith-Waterman traceback matrix.
/// Stores the I, D and S traceback matrix values in two bytes.
/// Possible traceback moves include : start, insert, delete, match, substitute,
/// prefix clip and suffix clip for x & y. So we need 4 bits each for matrices I, D, S
/// to keep track of these 9 moves.
#[derive(Copy, Clone)]
pub struct TracebackCell {
    v: u16,
}

impl Default for TracebackCell {
    fn default() -> Self {
        TracebackCell { v: 0 }
    }
}

// Traceback bit positions (LSB)
const I_POS: u8 = 0; // Meaning bits 0,1,2,3 corresponds to I and so on
const D_POS: u8 = 4;
const S_POS: u8 = 8;

// Traceback moves
const TB_START: u16 = 0b0000;
const TB_INS: u16 = 0b0001;
const TB_DEL: u16 = 0b0010;
const TB_SUBST: u16 = 0b0011;
const TB_MATCH: u16 = 0b0100;

const TB_XCLIP_PREFIX: u16 = 0b0101; // prefix clip of x
const TB_XCLIP_SUFFIX: u16 = 0b0110; // suffix clip of x
const TB_YCLIP_PREFIX: u16 = 0b0111; // prefix clip of y
const TB_YCLIP_SUFFIX: u16 = 0b1000; // suffix clip of y

const TB_MAX: u16 = 0b1000; // Useful in checking that the
                            // TB value we got is a valid one

impl TracebackCell {
    /// Initialize a blank traceback cell
    #[inline(always)]
    pub fn new() -> TracebackCell {
        Default::default()
    }

    /// Sets 4 bits [pos, pos+4) with the 4 LSBs of value
    #[inline(always)]
    fn set_bits(&mut self, pos: u8, value: u16) {
        let bits: u16 = (0b1111) << pos;
        assert!(
            value <= TB_MAX,
            "Expected a value <= TB_MAX while setting traceback bits"
        );
        self.v = (self.v & !bits) // First clear the bits
            | (value << pos) // And set the bits
    }

    #[inline(always)]
    pub fn set_i_bits(&mut self, value: u16) {
        // Traceback corresponding to matrix I
        self.set_bits(I_POS, value);
    }

    #[inline(always)]
    pub fn set_d_bits(&mut self, value: u16) {
        // Traceback corresponding to matrix D
        self.set_bits(D_POS, value);
    }

    #[inline(always)]
    pub fn set_s_bits(&mut self, value: u16) {
        // Traceback corresponding to matrix S
        self.set_bits(S_POS, value);
    }

    // Gets 4 bits [pos, pos+4) of v
    #[inline(always)]
    fn get_bits(self, pos: u8) -> u16 {
        (self.v >> pos) & (0b1111)
    }

    #[inline(always)]
    pub fn get_i_bits(self) -> u16 {
        self.get_bits(I_POS)
    }

    #[inline(always)]
    pub fn get_d_bits(self) -> u16 {
        self.get_bits(D_POS)
    }

    #[inline(always)]
    pub fn get_s_bits(self) -> u16 {
        self.get_bits(S_POS)
    }

    /// Set all matrices to the same value.
    pub fn set_all(&mut self, value: u16) {
        self.set_i_bits(value);
        self.set_d_bits(value);
        self.set_s_bits(value);
    }
}

/// Internal traceback.
struct Traceback {
    rows: usize,
    cols: usize,
    matrix: Vec<TracebackCell>,
}

impl Traceback {
    fn with_capacity(m: usize, n: usize) -> Self {
        let rows = m + 1;
        let cols = n + 1;
        Traceback {
            rows,
            cols,
            matrix: Vec::with_capacity(rows * cols),
        }
    }
    fn init(&mut self, m: usize, n: usize) {
        self.matrix.clear();
        let mut start = TracebackCell::new();
        start.set_all(TB_START);
        // set every cell to start
        self.resize(m, n, start);
    }

    #[inline(always)]
    fn set(&mut self, i: usize, j: usize, v: TracebackCell) {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        self.matrix[i * self.cols + j] = v;
    }

    #[inline(always)]
    fn get(&self, i: usize, j: usize) -> &TracebackCell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &self.matrix[i * self.cols + j]
    }

    fn get_mut(&mut self, i: usize, j: usize) -> &mut TracebackCell {
        debug_assert!(i < self.rows);
        debug_assert!(j < self.cols);
        &mut self.matrix[i * self.cols + j]
    }

    fn resize(&mut self, m: usize, n: usize, v: TracebackCell) {
        self.rows = m + 1;
        self.cols = n + 1;
        self.matrix.resize(self.rows * self.cols, v);
    }
}
