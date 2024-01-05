process MULTIQC {
    label 'process_medium'

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"

    input:
    path multiqc_config
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path methods_description
    path logo
    path ('fastqc/raw/*')
    path ('fastqc/trim/*')
    path ('fastp/*')
    path('alignment/rnalib/*')
    path('alignment/genome/*')
    path('alignment/nms/*')
    path('alignment/dedup/*')

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc \\
        -f \\
        $args \\
        $custom_config \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}