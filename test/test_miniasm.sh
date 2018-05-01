#!/usr/bin/env bash

nextflow run ../main.nf --cpus 2 --mem 4GB --reads ecoli_small.fasta --assembler miniasm --output test_miniasm
