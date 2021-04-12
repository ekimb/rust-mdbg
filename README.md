`rust-mdbg`: Minimizer-space *de Bruijn* graphs (mdBG) for whole-genome assembly
=========

`rust-mdbg` is an ultra-fast minimizer-space *de Bruijn* graph (mdBG) implementation, geared towards the assembly of long and accurate reads such as PacBio HiFi.

## Rationale

`rust-mdbg` performs mdBG construction of a 52x human genome HiFi data in under 30 minutes on 8 threads, with 16 GB of maximum RAM usage.

`rust-mdbg` is fast because it operates in *minimizer-space*, meaning that the reads, the assembly graph, and the final assembly, are all represented as ordered lists of minimizers, instead of strings of nucleotides. A conversion step then yields a classical *base-space* representation.

## Limitations

However, this high speed comes at a cost! :) 
* `rust-mdbg` gives good-quality results but still of lower contiguity and completeness than state-of-the-art assemblers such as HiCanu and hifiasm. 
* `rust-mdbg` performs best with at least 40x-50x of coverage.
* No polishing step is implemented; so, assemblies will have around the same accuracy as the reads.

## Installation

Clone the repository (make sure you have a working Rust environment), and run `cargo build --release`.

For performing graph simplifications, [gfatools](https://github.com/lh3/gfatools/) is required.

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
Usage:

`utils/multik <reads.fq.gz> <some_output_prefix> <nb_threads>`

## Overview

`rust-mdbg` is a **modular** assembler. It consists of three components:

 1) `rust-mdbg`, to perform assembly in minimizer-space
 2) `gfatools` (external component), to perform graph simplifications
 3) `to_basespace`, to convert a minimizer-space assembly to base-space

For convenience, components 2) and 3) are wrapped into a script called `magic_simplify`.

## Input

`rust-mdbg` takes a single FASTA/FASTQ input (gzip-compressed or not). Multi-line sequences are not supported.

For better results, please use the provided script [`utils/remove_homopoly.py`](https://raw.githubusercontent.com/ekimb/rust-mdbg/master/utils/remove_homopoly.py) to homopolymer-compress your reads. 

Usage example:

`seqtk seq -A reads.fq.gz | python remove_homopoly.py /dev/stdin | gzip -1 -c > reads.fa.hpc.gz`

## Output data 

The output of `rust-mdbg` consists of:

* A `.gfa` file containing the minimizer-space *de Bruijn* graph, without sequences,
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

The main parameters of `rust-mdbg` are the `k`-min-mer value `k`, the minimizer length `l`, and the minimizer density `d` (delta in the paper). 
For better results, and also without the need to set any parameter, try the multi-`k` strategy (see Multi-`k` assembly section). 
This section explains how parameters are set in single-`k` assembly.

All three parameters `k`, `l`, and `d` significantly impact the quality of results. One can think of them as a generalization of the `k` parameter in classical *de Bruijn* graphs. When you run `rust-mdbg` without specifying parameters, it sets them to:

   `d` = 0.003

   `l` = 12

   `k` = 0.75 * `average_readlen` * `d`
   
These parameters will give reasonable, but far from optimal, draft assemblies. We experimentally found that the best results are often obtained with `k` values within 20-40, `l` within 10-14, and `d` within 0.001-0.005. Setting `k` and `d` such that the ratio `k`/`d` is slightly below the read length appears to be an effective strategy. 


## Performance

|Dataset                 | Genome size (HPC)   | Coverage  | <div style="width:1090px">Parameters</div> | N50     | Runtime | Memory |
|:-----------------------|:-------------:|:----:|------------------------------------:|--------:|:------------------------------------------|-------:|
|[*D. melanogaster* HiFi](http://www.ncbi.nlm.nih.gov/bioproject/?term=SRR10238607)    | 98 Mbp | 100x | auto<br>multi-k<br>k=35,l=12,d=0.002 | 2.5Mbp<br>2.5Mbp<br>3.9Mbp  |  2m15s<br>15mins<br>1m40s                  |   2.5G<br>1.8G<br>6.7G |
|[Strawberry HiFi](http://www.ncbi.nlm.nih.gov/bioproject/?term=SRR11606867)    | 0.7 Gbp | 36x | auto<br>multi-k<br>k=38,l=14,d=0.003| 0.5Mbp<br>1Mbp<br>0.7Mbp  |  6m12s<br>40mins<br>5m31s                  |   12G<br>11G<br>10G |
|[*H. sapiens* (HG002) HiFi](https://github.com/human-pangenomics/HG002_Data_Freeze_v1.0#pacbio-hifi-1)  | 2.2 Gbp | 52x  | auto<br>multi-`k`<br>k=21,l=14,d=0.003 | 1.0Mbp<br>16.9Mbp<br>13.6Mbp |  27m30s<br>3h15m<br>24m47s           | 16.9G<br>20G<br>10.6G |

Runtime breakdown:<br>
*D. melanogaster*: 1m40s = 1m18s  `rust-mdbg` + 8s `gfatools` + 14s  `to_basespace`<br>
*H. sapiens*: 24m47s = 18m58s `rust-mdbg` + 3m19s `gfatools` + 2m30s `to_basespace`

The runs with custom parameters (from the paper) were made with commit `b99d938`, and unlike in the paper, we did not use robust minimizers which requires additional `l`-mer counting beforehand.
Reads and assemblies were homopolymer-compressed and the homopolymer-compressed genome size is reported.
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

* [Barış Ekim](http://people.csail.mit.edu/ekim/), supervised by [Bonnie Berger](http://people.csail.mit.edu/bab/) at the Computer Science and Artificial Intelligence Laboratory (CSAIL) at Massachusetts Institute of Technology (MIT), and
* [Rayan Chikhi](http://rayan.chikhi.name) at the Department of Computational Biology at Institut Pasteur.

## Citation

A pre-print will be posted at a later date.

## Contact

Should you have any inquiries, please contact [Barış Ekim](http://people.csail.mit.edu/ekim/) at baris [at] mit [dot] edu, or [Rayan Chikhi](http://rayan.chikhi.name) at rchikhi [at] pasteur [dot] fr.


