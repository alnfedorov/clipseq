process ANNOTATE_XLINKS {
    tag "$meta.id"

    // conda "conda-forge::python=3.10.11 bioconda::pysam=0.21.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py310h41dec4a_1' :
    //     'biocontainers/pysam:0.21.0--py310h41dec4a_1' }"

    input:
    tuple val(meta), path(bam, stageAs: "input.bam"), path(bai, stageAs: "input.bam.bai")
    path(rnalib_manifest, stageAs: "rnalib_manifest.csv.gz")
    path(genomic_regions, stageAs: "genomic_regions.bed.gz")
    val unclassified_label
    val summary_level

    output:
    path "${meta.id}.counts.csv.gz"           , emit: counts
    path "${meta.id}.counts-summary.tsv"      , emit: counts_summary
    path "${meta.id}.collisions.csv.gz"       , emit: collisions
    path "${meta.id}.collisions-summary.tsv"  , emit: collisions_summary

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    annotate_xlink.py \\
        --threads ${task.cpus} \\
        --bam ${bam} \\
        --rnalib-manifest ${rnalib_manifest} \\
        --genomic-regions ${genomic_regions} \\
        --unclassified-label ${unclassified_label} \\
        --summary-level ${summary_level} \\
        --output-counts ${counts} \\
        --output-counts-summary ${counts_summary} \\
        --output-collisions ${collisions} \\
        --output-collisions-summary ${collisions_summary}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: \$(python -c 'import pysam; print(pysam.__version__)')
        intervaltree: \$(python -c 'import intervaltree; print(intervaltree.__version__)')
        joblib: \$(python -c 'import joblib; print(joblib.__version__)')
    END_VERSIONS
    """
}
