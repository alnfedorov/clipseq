process NMS {
    tag "$meta.id"

    conda "conda-forge::python=3.10.11 bioconda::pysam=0.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py310h41dec4a_1' :
        'biocontainers/pysam:0.21.0--py310h41dec4a_1' }"

    input:
    tuple val(meta), path(reads, stageAs: "reads.fastq.gz"), path(rnalib, stageAs: "rnalib.bam"), path(reads, stageAs: "genome.bam")
    
    output:
    tuple val(meta), path("${meta.id}.nms.bam"), path("${meta.id}.nms.bam.bai"), emit: alignment 
    path "${meta.id}.unmapped.fastq.gz"                                        , emit: unmapped
    path "${meta.id}_nms.tsv"                                                  , emit: summary
    path "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    nms.py \\
        ${task.cpus} \\
        reads.fastq.gz \\
        rnalib.bam \\
        genome.bam \\
        ${meta.id}_nms.tsv \\
        ${meta.id}.nms.bam \\
        ${meta.id}.unmapped.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
}
