# Nanoflow

[![build-status](https://img.shields.io/travis/HadrienG/nanoflow/master.svg?style=flat-square)](https://travis-ci.org/HadrienG/nanoflow)
[![made-with-nextflow](https://img.shields.io/badge/made%20with-nextflow-green.svg?longCache=true&style=flat-square)](https://www.nextflow.io/)
[![gpl](https://img.shields.io/badge/license-GPL--3.0-lightgrey.svg?style=flat-square)](LICENSE)

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

You'll also need docker installed if you wish to run the pipeline locally.

## Usage

Once you have nextflow and docker installed:

```bash
nextflow run hadrieng/nanoflow --reads reads.fasta --assembler unicycler --output results
```

### Options

#### --cpus
* number of cpus
* 4 by default

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

### profiles

By default nanoflow will execute on your local nachine using dockers.
You can modify this behavior with the `-profile` option

Existing profiles:

* `planet`: to execute it on the SGBC cluster

## License

Code is [GPL-3.0](LICENSE)

## Contributing

We welcome contributions from the community! See our
[Contributing](CONTRIBUTING.md) guidelines
