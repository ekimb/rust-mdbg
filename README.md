`rust-mdbg`: Minimizer-space de Bruijn graphs (mdBG) for whole-genome assembly
=========

`rust-mdbg` is an ultra-fast minimizer-space de Bruijn graph (mdBG) implementation, geared towards the assembly of long and accurate reads such as PacBio HiFi.

## Overview and purpose

`rust-mdbg` does not replace a conventional genome assembler, as it has slightly lower contiguity and completeness. However, it can be useful
for quickly assembling reads. It performs mdBG construction of a human genome in 26 minutes on 8 threads, with 16 GB of maximum RAM usage.

`rust-mdbg` operates in *minimizer-space*, meaning that the reads, the assembly graph, and the final assembly, are all seen as ordered lists of minimizers, instead of strings of nucleotides. A conversion step is necessary to obtain a classical *base-space* representation.


## Installation

Clone the repository (make sure you have a working Rust environment), and run `cargo build --release`.

For performing graph simplifications, [gfatools](https://github.com/lh3/gfatools/) is required.

## Quick start

```
cargo build --release
target/release/rust-mdbg reads-0.00.fa.gz -k 7 --density 0.0008 -l 10 --minabund 2 --prefix example
utils/magic_simplify example
```

## Overview

`rust-mdbg` is a **modular** assembler. It consists of three components:

 1) `rust-mdbg`, to perform assembly in minimizer-space
 2) `gfatools` (external component), to perform graph simplifications
 3) `to_basespace`, to convert a minimizer-space assembly to base-space

For convenience, components 2) and 3) are wrapped into a script called `magic_simplify`.

## Input

`rust-mdbg` takes a single FASTA/FASTQ input (gzip-compressed or not), and requires the user to specify the k, density and l values, as discussed in the article. 

## Output data

The output of the `rust-mdbg` program consists of:

* a `.gfa` file containing the minimizer-space de Bruijn graph, without sequences
* a `.sequences` file containing the sequences of the nodes of the graph

The `to_basespace` program allows to combine both outputs and produde a `.gfa` file with sequences.

## Running an example

A sample set of reads was provided in the `example/` folder. Run

`target/release/rust-mdbg reads-0.00.fa.gz -k 7 --density 0.0008 -l 10 --minabund 2 --prefix example`

which will create an `example.gfa` file.

In order to populate the `.gfa` file with base-space sequences and perform graph simplification, run

`utils/magic_simplify example`

which will create `example.msimpl.gfa` and `example.msimpl.fa` files.

In the case that you only want to convert to base-space with no graph simplification, there are two ways:
* with `gfatools`

```
gfatools asm -u  example.gfa > example.unitigs.gfa
target/release/to_basespace --gfa example.unitigs.gfa --sequences example.sequences
```

* without `gfatools` (slower, but the code is more straightforward to understand)

`utils/complete_gfa.py example.sequences example.gfa`

In both cases this will create an `example.complete.gfa` file that you can convert to FASTA with

`bash $DIR/gfa2fasta.sh example.complete`

## Parameters

The main parameters of `rust-mdbg` are the k-min-mer value `k`, the minimizer length `l`, and the minimizer density `d` (delta in the paper).
Have a look at the next section for some examples of how to set them. For optimal results, we recommend trying `k` values within 20-40, `l` within 10-14, and `d` within 0.001-0.005.

## Performance

|Dataset                 | Genome size   | Cov  | Parameters                           | N50     | Time (rust-mdbg + gfatools + to_basespace) | Memory |
|:-----------------------|:-------------:|:----:|:------------------------------------:|--------:|:------------------------------------------|-------:|
|[D. melanogaster HiFi](http://www.ncbi.nlm.nih.gov/bioproject/?term=SRR10238607)    | 130 Mbp | 100x | k=35 l=12 d=0.002  | 3.9Mbp  |  1m40s (1m18s + 8s + 14s)                  |   1.8G |
|[H. Sapiens HG002 HiFi Sequel II chem 2.0](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0#pacbio-hifi-1)  | 2.2 Gbp | 52x  | k=21 l=14 d=0.003 presimp=0.01 | 13.6Mbp |  21m48s (15m59s + 3m19s + 2m30s)           |  14.5G |

The following runs were made with commit `b99d938`, and unlike in the paper we did not use robust minimizers which requires additional l-mer counting beforehand.
Reads were homopolymer-compressed and the genome size is also the homopolymer-compressed one.
In addition to the parameters shown in the table, the `rust-mdbg` command line also contained: `--bf --no-error-correct --threads 8`.

## License

`rust-mdbg` is freely available under the [MIT License](https://opensource.org/licenses/MIT).

## Contributors

Development of `rust-mdbg` is led by [Barış Ekim](http://people.csail.mit.edu/ekim/), collaboratively in the labs of [Bonnie Berger](http://people.csail.mit.edu/bab/) at the Computer Science and Artificial Intelligence Laboratory (CSAIL) at Massachusetts Institute of Technology (MIT), and [Rayan Chikhi](http://rayan.chikhi.name) at the Department of Computational Biology at Institut Pasteur.

## Citation

A pre-print will be posted at a later date.

## Contact

Should you have any inquiries, please contact [Barış Ekim](http://people.csail.mit.edu/ekim/) at baris [at] mit [dot] edu, or [Rayan Chikhi](http://rayan.chikhi.name) at rchikhi [at] pasteur [dot] fr.


