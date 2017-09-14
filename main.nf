#!/usr/bin/env nextflow

params.reads = ''
params.assembler = ''
params.consensus = ''

process adapter_trimming {
	
    input:
	file(reads) from Channel.value( file(params.reads) )
	
    output:
	file("trimmed.fastq") into trimmedreads
	
	script:
    """	
	echo ${task.cpus} ${task.memory} ${task.jobName}
	#Adapters are trimmed with porechop
	porechop -i $reads -t ${task.cpus} -o trimmed.fastq	
    """
}

process assembly {
    input:

    output:

    """
    """
}

process consensus {
    input:

    output:

    """
    """
}
