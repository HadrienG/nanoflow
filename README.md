# Nanoflow

[![build-status](https://img.shields.io/travis/HadrienG/nanoflow/master.svg?style=flat-square)](https://travis-ci.org/HadrienG/nanoflow)
[![made-with-nextflow](https://img.shields.io/badge/made%20with-nextflow-green.svg?longCache=true&style=flat-square)](https://www.nextflow.io/)
[![gpl](https://img.shields.io/badge/license-GPL--3.0-lightgrey.svg?style=flat-square)](LICENSE)

## Introduction

Nanoflow aims at being a reproducible pipeline for bacterial genome assembly
of nanopore reads.

Nanoflow uses the following software

| Task | Software | Version | Docker | Lmod + SGE |
| --- | --- | --- | --- | --- |
| adapter trimming | [porechop](https://github.com/rrwick/Porechop) | 0.2.3 | ![Docker Build Status](https://img.shields.io/docker/build/hadrieng/porechop.svg?style=flat-square) | ✔️ |
| assembly | [miniasm](https://github.com/lh3/miniasm) + [minimap2](https://github.com/lh3/minimap2) | 0.2-r168 / 2.10-r768 | ![Docker Build Status](https://img.shields.io/docker/build/hadrieng/miniasm.svg?style=flat-square) | ✔️ |
| | [canu](https://github.com/marbl/canu) | Unknown | ![Docker Build Status](https://img.shields.io/docker/build/hadrieng/canu.svg?style=flat-square) | ✔️ |
| | [unicycler](https://github.com/rrwick/Unicycler) | 0.4.5 | ![Docker Build Status](https://img.shields.io/docker/build/hadrieng/unicycler.svg?style=flat-square) | ✔️ |
| consensus | [racon](https://github.com/isovic/racon) + [minimap2](https://github.com/lh3/minimap2) | 1.3.1 / 2.10-r768 |  ![Docker Build Status](https://img.shields.io/docker/build/hadrieng/racon.svg?style=flat-square) | ✔️ |
| polishing | [nanopolish](https://github.com/jts/nanopolish) | Unknown | ![Docker Build Status](https://img.shields.io/docker/build/hadrieng/nanopolish.svg?style=flat-square) | ✔️ |

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

#### --mem
* amount of memory to use, in GB
* 8GB by default
* has to end with "GB"

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
