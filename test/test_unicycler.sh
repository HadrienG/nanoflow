#!/usr/bin/env bash

nextflow run ../main.nf --cpus 2 --reads ecoli_small.fasta --assembler unicycler --output test_results
