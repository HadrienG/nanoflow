#!/usr/bin/env nextflow

params.reads = ''
params.genome_size = ''

params.workingdir = 'results'

// choose the assembler
params.assembler = 'miniasm'
if (params.assembler != 'miniasm' && params.assembler != 'canu') {
    exit 1, "Error --assembler: ${params.assembler}. \
    Should be 'miniasm' or 'canu'"
}
// requires genome_size for canu
if (params.assembler == 'canu' && params.genome_size == '') {
    exit 1, "Error --genome_size is a required paramater for canu"
}

process adapter_trimming {
    input:
	file(reads) from Channel.value( file(params.reads) )

    output:
	file("trimmed.fastq") into trimmed_reads

	script:
    """
	echo ${task.cpus} ${task.memory}
	#Adapters are trimmed with porechop
	porechop -i $reads -t ${task.cpus} -o trimmed.fastq
    """
}

// Trimmed reads are used by both 'assembly' and 'consensus',
// channel trimmed_reads is duplicated
trimmed_reads.into { trimmed_for_assembly; trimmed_for_consensus }

process assembly {
    publishDir params.workingdir, mode: 'copy', pattern: "assembly.fasta"

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
	publishDir params.workingdir, mode: 'copy', pattern: "assembly_consensus.fasta"

    input:
	file(reads) from trimmed_for_consensus
	file(assembly) from assembly

    output:
	file 'assembly_consensus.fasta' into assembly_consensus

	script:
	if(params.assembler == 'miniasm')
    """
	minimap2 -x map-ont -t ${task.cpus} ${assembly} ${reads} > assembly.paf
	racon -t ${task.cpus} ${reads} assembly.paf ${assembly} assembly_consensus.fasta
	"""
    else if(params.assembler == 'canu')
    """
    echo \"Consensus was already performed by canu\"
    cp ${assembly} assembly_consensus.fasta
    """
}
