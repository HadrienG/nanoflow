#!/usr/bin/env nextflow

params.fast5 = ''
params.reads = ''
params.genome_size = ''

params.output = 'results'

// choose the assembler
params.assembler = 'miniasm'
if (params.assembler != 'miniasm' && params.assembler != 'canu') {
    exit 1, "Error --assembler: ${params.assembler}. \
    Should be 'miniasm' or 'canu'"
}
// requires genome_size for canu
if (params.assembler == 'canu' && params.genome_size == '') {
    exit 1, 'Error --genome_size is a required paramater for canu'
}

process adapter_trimming {
    input:
	   file(reads) from file(params.reads)

    output:
	   file('trimmed.fastq') into trimmed_reads

	script:
        """
    	porechop -i "${reads}" -t "${task.cpus}" -o trimmed.fastq
        """
}

// Trimmed reads are used by both 'assembly' and 'consensus',
// channel trimmed_reads is duplicated
trimmed_reads.into { trimmed_for_assembly; trimmed_for_consensus }

process assembly {
    publishDir params.output, mode: 'copy', pattern: 'assembly.fasta'

    input:
        file reads from trimmed_for_assembly
        val genome_size from params.genome_size
    output:
        file 'assembly.fasta' into assembly

    script:
    if(params.assembler == 'miniasm')
        """
        minimap2 -x ava-ont -t "${task.cpus}" "${reads}" "${reads}" > "${reads}.paf"
        miniasm -f "${reads}" "${reads}.paf" > "${reads}.gfa"
        awk '/^S/{print ">"\$2"\\n"\$3}' "${reads}.gfa" | fold > assembly.fasta
        """
    else if(params.assembler == 'canu')
        """
        canu -p assembly -d canu_out \
            genomeSize="${genome_size}" -nanopore-raw "${reads}" \
            maxThreads="${task.cpus}" useGrid=false gnuplotTested=true
        mv canu_out/assembly.contigs.fasta assembly.fasta
        """
}


process consensus {
	publishDir params.output, mode: 'copy', pattern: 'assembly_consensus.fasta'
    echo true

    input:
    	file(reads) from trimmed_for_consensus
    	file(assembly) from assembly

    output:
	   file 'assembly_consensus.fasta' into assembly_consensus

	script:
	if(params.assembler == 'miniasm')
        """
    	minimap2 -x map-ont -t "${task.cpus}" "${assembly}" "${reads}" > assembly.paf
    	racon -t "${task.cpus}" "${reads}" assembly.paf "${assembly}" assembly_consensus.fasta
    	"""
    else if(params.assembler == 'canu')
        """
        echo "[info] consensus was already performed by canu. Skipping."
        cp "${assembly}" assembly_consensus.fasta
        """
}


process polishing {
    publishDir params.output, mode: 'copy', pattern: 'polished_genome.fa'
    echo true

    input:
        file(assembly) from assembly_consensus
        file(reads) from file(params.reads)
        val(fast5_dir) from params.fast5

    output:
        file 'polished_genome.fa' into assembly_polished

    script:
    if (params.fast5 != '')
        """
        nanopolish index -d "${fast5_dir}" "${reads}"
        bwa index "${assembly}"
        bwa mem -x ont2d -t "${task.cpus}" "${assembly}" "${reads}" | \
            samtools sort -o reads.sorted.bam -T reads.tmp -
        samtools index reads.sorted.bam
        nanopolish_makerange.py "${assembly}" | parallel --results \
            nanopolish.results -P "${task.cpus}" nanopolish variants --consensus \
            polished.{1}.fa -w {1} -r "${reads}" -b reads.sorted.bam -g \
            "${assembly}" -t 1 --min-candidate-frequency 0.1
        nanopolish_merge.py polished.*.fa > polished_genome.fa
        """
    else
        """
        echo "[warn] polishing requires fast5 data. Skipping."
        cp "${assembly}" polished_genome.fa
        """
}
