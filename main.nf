#!/usr/bin/env nextflow

params.reads = ''
params.genome_size = ''

// choose the assembler
params.assembler = 'miniasm'
if (params.assembler != 'miniasm' && params.assembler != 'canu'){
    exit 1, "Error --assembler: ${params.assembler}. \
    Should be 'miniasm' or 'canu'"
}
params.consensus = ''

process adapter_trimming {

    input:
	file(reads) from Channel.value( file(params.reads) )

    output:
	file("trimmed.fastq") into trimmed_reads

	script:
    """
	echo ${task.cpus} ${task.memory} ${task.jobName}
	#Adapters are trimmed with porechop
	porechop -i $reads -t ${task.cpus} -o trimmed.fastq
    """
}

process assembly {
    publishDir 'results'
    container {params.assembler == 'miniasm' ? 'hadrieng/miniasm' : 'hadrieng/canu'}
    cpus 2

    input:
        file reads from trimmed_reads
        val genome_size from params.genome_size
    output:
        file 'assembly.fasta' into assembly

    script:
    if(params.assembler == 'miniasm')
        """
        minimap2 -x ava-ont -t "${task.cpus}" "${reads}" "${reads}" | \
        gzip -1 > "${reads}.paf.gz"
        miniasm -f "${reads}" "${reads}.paf.gz" > \
        "${reads}.gfa"
        awk '/^S/{print ">"\$2"\\n"\$3}' "${reads}.gfa" | \
        fold > assembly.fasta
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
    input:

    output:

    """
    """
}
