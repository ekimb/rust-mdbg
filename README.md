rust-mdbg
=========

`rust-mdbg` is an ultra-fast minimizer-space de Bruijn graph (mdBG) implementation, geared towards the assembly of long and accurate reads such as PacBio HiFi.

# Purpose

`rust-mdbg` does not replace a conventional genome assembler, as it has slightly lower contiguity and completeness. However, it can be useful
for quickly assembling reads. It performs mDBG construction of a human genome in 22 minutes on 8 threads, with 9 GB of maximum RAM usage.

# Installation

Clone the repository (make sure you have a working Rust environment), and run `cargo build --release`.

For performing graph simplifications, [gfatools](https://github.com/lh3/gfatools/) is required.

# Quick start

Currently, `rust-mdbg` only takes a single FASTA input, and requires the user to specify the k, density and l values, as discussed in the article. 

A sample set of reads was provided in the `example/` folder.  Uncompress them:

`gunzip example/reads-0.00.fa.gz`

then run:

`target/release/rust-mdbg reads-0.00.fa -k 7 --density 0.0008 -l 10 --minabund 2 --prefix`

which will create an `example.gfa` file.

In order to populate the `.gfa` file with base-space sequences and perform graph simplification, run

`utils/magic_simplify example`

which will create `example.msimpl.gfa` and `example.msimpl.fa` files.

In the case that you only want to convert to base-space with no graph simplification, run

`target/release/to_basespace --gfa example.gfa --sequences example.sequences`

which will create an `example.complete.gfa` file that you can convert to FASTA with

`bash $DIR/gfa2fasta.sh example.complete`
