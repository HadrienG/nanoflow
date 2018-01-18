# Nanoflow

## Introduction

Nanoflow aims at being a reproducible pipeline for bacterial genome assembly
of nanopore reads.

Nanoflow is a work in progress. At the moment it offers:

- [x] adapter trimming using [porechop](https://github.com/rrwick/Porechop)
- [x] assembly using [miniasm](https://github.com/lh3/miniasm), [canu](https://github.com/marbl/canu) or [unicycler](https://github.com/rrwick/Unicycler)
- [x] consensus using [racon](https://github.com/isovic/racon)
- [x] polishing using [nanopolish](https://github.com/jts/nanopolish)

## Installation

To install nextflow (make sure you have java installed):

```bash
curl -fsSL get.nextflow.io | bash
```

You'll also need docker installed.

## Usage

Once you have nextflow and docker installed:

```bash
nextflow run hadrieng/nanoflow --reads reads.fasta --assembler unicycler --output results
```

### Options

#### --assembler
* which assembler to use
* can be miniasm, unicycler or canu
* required

#### --reads
* location of the nanopore reads in fasta or fastq
* required

#### --output
* output directory
* required

#### --fast5
* location of the fast5 files
* required for polishing with nanopolish

#### --genome_size
* genome size of the organism you are trying to assemble
* required if `--assembler canu`

## License

Code is [GPL-3.0](LICENSE)

## Contributing

We welcome contributions from the community! See our
[Contributing](CONTRIBUTING.md) guidelines
