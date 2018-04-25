#!/usr/bin/env bash

nextflow run ../main.nf --reads ecoli_small.fasta --assembler unicycler --output test_results
