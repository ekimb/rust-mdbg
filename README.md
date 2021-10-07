# `rust-mdbg`: Minimizer-space de Bruijn graphs (mdBG) for whole-genome assembly <br>
[![DOI](https://zenodo.org/badge/310619686.svg)](https://zenodo.org/badge/latestdoi/310619686)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/ekimb/rust-mdbg)
![GitHub last commit](https://img.shields.io/github/last-commit/ekimb/rust-mdbg)
![GitHub](https://img.shields.io/github/license/ekimb/rust-mdbg)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/rust-mdbg/README.html)

`rust-mdbg` is an ultra-fast minimizer-space de Bruijn graph (mdBG) implementation, geared towards the assembly of long and accurate reads such as [PacBio HiFi](https://www.pacb.com/smrt-science/smrt-sequencing/hifi-reads-for-highly-accurate-long-read-sequencing/).

## Rationale

`rust-mdbg` performs mdBG construction of a 52x human genome HiFi data in around 10 minutes on 8 threads, with 10GB of maximum RAM usage.

`rust-mdbg` is fast because it operates in *minimizer-space*, meaning that the reads, the assembly graph, and the final assembly, are all represented as ordered lists of minimizers, instead of strings of nucleotides. A conversion step then yields a classical *base-space* representation.

## Limitations

However, this high speed comes at a cost! :) 
* `rust-mdbg` gives good-quality results but still of lower contiguity and completeness than state-of-the-art assemblers such as [`HiCanu`](https://github.com/marbl/canu) and [`hifiasm`](https://github.com/chhylp123/hifiasm). 
* `rust-mdbg` performs best with at least 40x to 50x of coverage.
* No polishing step is implemented; so, assemblies will have around the same accuracy as the reads.
* Cannot assemble Nanopore data due to its higher error rate (see [this comment](https://github.com/ekimb/rust-mdbg/issues/4#issuecomment-860817828))

## Installation

Clone the repository (make sure you have a working Rust environment), and run 

`cargo build --release`

Alternatively, you can install from [`bioconda`](https://bioconda.github.io/index.html):

`conda install -c bioconda rust-mdbg`

which has the Rust binaries, but not the additional scripts. For performing graph simplifications, [`gfatools`](https://github.com/lh3/gfatools/) is required.

## Quick start

```
cargo build --release
target/release/rust-mdbg reads-0.00.fa.gz -k 7 --density 0.0008 -l 10 --minabund 2 --prefix example
utils/magic_simplify example
```

## Multi-`k` assembly

For better contiguity, try the provided multi-`k` assembly script.
It performs assembly iteratively, starting with `k`= 10, up to an automatically-determined largest `k`. 
This comes at the expense of ~7x longer running time.

`utils/multik <reads.fq.gz> <some_output_prefix> <nb_threads>`

## Overview

`rust-mdbg` is a **modular** assembler. It consists of three components:

 1) `rust-mdbg`, to perform assembly in minimizer-space
 2) `gfatools` (external component), to perform graph simplifications
 3) `to_basespace`, to convert a minimizer-space assembly to base-space

For convenience, components 2 and 3 are wrapped into a script called `magic_simplify`.

## Input

`rust-mdbg` takes a single FASTA/FASTQ input (gzip-compressed or not). Multi-line sequences, and sequences with lowercase characters, are not supported. 

If you have [`seqtk`](https://github.com/lh3/seqtk) installed, you can use

`seqtk seq -AU reads.unformatted.fq > reads.fa`

to format reads accordingly.

## Output data 

The output of `rust-mdbg` consists of:

* A `.gfa` file containing the minimizer-space de Bruijn graph, without sequences,
* Several `.sequences` files containing the sequences of the nodes of the graph.

The executable `to_basespace` allows to combine both outputs and produce a `.gfa` file, with sequences.

## Running an example

A sample set of reads is provided in the `example/` folder. Run

`target/release/rust-mdbg reads-0.00.fa.gz -k 7 --density 0.0008 -l 10 --minabund 2 --prefix example`

which will create an `example.gfa` file.

In order to populate the `.gfa` file with base-space sequences and perform graph simplification, run

`utils/magic_simplify example`

which will create `example.msimpl.gfa` and `example.msimpl.fa` files.


## Parameters

The main parameters of `rust-mdbg` are the `k`-min-mer value `k`, the minimizer length `l`, and the minimizer density `d` (delta in the paper).  Another parameter is `--presimp`, set by default to 0.01, which performs a graph simplification: a neighbor node is deleted if its abundance is below 1% that of `min(max(abundance of other neighbors), abundance of current node)`.
For better results, and also without the need to set any parameter, try the multi-`k` strategy (see Multi-`k` assembly section). 
This section explains how parameters are set in single-`k` assembly.

All three parameters `k`, `l`, and `d` significantly impact the quality of results. One can think of them as a generalization of the `k` parameter in classical de Bruijn graphs. When you run `rust-mdbg` without specifying parameters, it sets them to:

   `d` = 0.003

   `l` = 12

   `k` = 0.75 * `average_readlen` * `d`
   
These parameters will give reasonable, but far from optimal, draft assemblies. We experimentally found that the best results are often obtained with `k` values within 20-40, `l` within 10-14, and `d` within 0.001-0.005. Setting `k` and `d` such that the ratio `k`/`d` is slightly below the read length appears to be an effective strategy.

For further information on usage and parameters, run

`target/release/rust-mdbg -h`

for a one-line summary of each flag, or run

`target/release/rust-mdbg --help`

for a lengthy explanation of each flag.

## Performance

|Dataset                 | Genome size (HPC)   | Coverage  | <div style="width:1200px">Parameters</div> | N50     | Runtime | Memory |
|:-----------------------|:-------------:|:----:|------------------------------------:|--------:|:------------------------------------------|-------:|
|[*D. melanogaster* HiFi](http://www.ncbi.nlm.nih.gov/bioproject/?term=SRR10238607)    | 98Mbp | 100x | auto<br>multi-`k`<br>`k`=35,`l`=12,`d`=0.002 | 2.5Mbp<br>2.5Mbp<br>6.0Mbp  |  2m15s<br>15m<br>1m9s                  |   2.5GB<br>1.8GB<br>1.5GB |
|[Strawberry HiFi](http://www.ncbi.nlm.nih.gov/bioproject/?term=SRR11606867)    | 0.7Gbp | 36x | auto<br>multi-`k`<br>`k`=38,`l`=14,`d`=0.003| 0.5Mbp<br>1Mbp<br>0.7Mbp  |  6m12s<br>40m<br>5m31s                  |   12GB<br>11GB<br>10GB |
|[*H. sapiens* (HG002) HiFi](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0#pacbio-hifi-1)  | 2.2Gbp | 52x  | auto<br>multi-`k`<br>`k`=21,`l`=14,`d`=0.003 | 1.0Mbp<br>16.9Mbp<br>13.9Mbp |  27m30s<br>3h15m<br>10m23s           | 16.9GB<br>20GB<br>10.1GB |

Runtime breakdown:<br>
*H. sapiens*: 10m23s = 6m51s `rust-mdbg` + 1m48s `gfatools` + 1m44s `to_basespace`

The runs with custom parameters (from the paper) were made with commit `b99d938`, and unlike in the paper, we did not use robust minimizers which requires additional `l`-mer counting beforehand.
For historical reasons, reads and assemblies were homopolymer-compressed in those experiments and the homopolymer-compressed genome size is reported. So beware that these numbers are not directly comparable to the output of other assemblers.
In addition to the parameters shown in the table, the `rust-mdbg` command line also contained `--bf --no-error-correct --threads 8`.

## Running `rust-mdbg` without graph simplifications

To convert an assembly to base-space without performing any graph simplifications, there are two ways:

* with `gfatools`

```
gfatools asm -u  example.gfa > example.unitigs.gfa
target/release/to_basespace --gfa example.unitigs.gfa --sequences example.sequences
```

* without `gfatools` (slower, but the code is more straightforward to understand)

`utils/complete_gfa.py example.sequences example.gfa`

In both cases, this will create an `example.complete.gfa` file that you can convert to FASTA with

`bash utils/gfa2fasta.sh example.complete`

## License

`rust-mdbg` is freely available under the [MIT License](https://opensource.org/licenses/MIT).

## Developers

* [Barış Ekim](http://people.csail.mit.edu/ekim/), supervised by [Bonnie Berger](http://people.csail.mit.edu/bab/) at the Computer Science and Artificial Intelligence Laboratory (CSAIL) at Massachusetts Institute of Technology (MIT)
* [Rayan Chikhi](http://rayan.chikhi.name) at the Department of Computational Biology at Institut Pasteur

## Citation
* Barış Ekim, Bonnie Berger, and Rayan Chikhi, [Minimizer-space de Bruijn graphs: Whole-genome assembly of long reads in minutes on a personal computer](https://www.sciencedirect.com/science/article/pii/S240547122100332X), Cell Systems (2021).
* Barış Ekim, Bonnie Berger, and Rayan Chikhi, [Minimizer-space de Bruijn graphs](https://www.biorxiv.org/content/10.1101/2021.06.09.447586v1), biorXiv preprint (2021).
```
@article {mdbg,
	author = {Ekim, Bar{\i}{\c s} and Berger, Bonnie and Chikhi, Rayan},
	title = {Minimizer-space de Bruijn graphs: Whole-genome assembly of long reads in minutes on a personal computer},
	journal = {Cell Systems},
	year = {2021},
	issn = {2405-4712},
	doi = {https://doi.org/10.1016/j.cels.2021.08.009}
}
```


## Contact

Should you have any inquiries, please contact [Barış Ekim](http://people.csail.mit.edu/ekim/) at baris [at] mit [dot] edu, or [Rayan Chikhi](http://rayan.chikhi.name) at rchikhi [at] pasteur [dot] fr.


