process FIX_BAM_FLAGS_TAGS {
    tag "$meta.id"

    conda "conda-forge::python=3.10.11 bioconda::pysam=0.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py310h41dec4a_1' :
        'biocontainers/pysam:0.21.0--py310h41dec4a_1' }"

    input:
    tuple val(meta), path(bam, stageAs: "input.bam"), path(bai, stageAs: "input.bam.bai")
    
    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: alignment 
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir output

    fix-flags.py ${task.cpus} input.bam ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
    END_VERSIONS
    """
}
