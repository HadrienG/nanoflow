#!/usr/bin/env nextflow

params.reads = ''
params.workingdir = 'results'

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

// Trimmed reads are used by both 'assembly' and 'consensus', 
// channel trimmed_reads is duplicated
trimmed_reads.into { trimmed_for_assembly; trimmed_for_consensus }

process assembly {
    publishDir 'results'
    container 'hadrieng/miniasm'
    cpus 2

    input:
        file reads from trimmed_for_assembly
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
        echo TODO
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
	minimap -x map10k -t ${task.cpus} ${assembly} ${reads} assembly.paf
	racon -t ${task.cpus} ${reads} assembly.paf $assembly assembly_consensus.fasta	
	"""
    else if(params.assembler == 'canu')
    """
    echo \"Consensus was already performed by canu\"
    """
}

process annotation {
       container 'hadrieng/prokka'
       publishDir 'results'

    input:
        file input from assembly_consensus

    output:
        file "prokka/${input}.gff" into annotation_gff

    """
    prokka --outdir prokka --prefix ${input.baseName} --kingdom Bacteria $input
    """
}
