process MAKE_XLINK_TRACKS {
    tag "$meta.id"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f01e242bdea19948f0576fdca94777242fe4c2cb:a2983f328c2d85ae07c18ddf6651b25cf7d55187-0' :
        'biocontainers/mulled-v2-f01e242bdea19948f0576fdca94777242fe4c2cb:a2983f328c2d85ae07c18ddf6651b25cf7d55187-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.xln.bed.gz")       , emit: bed
    tuple val(meta), path("${meta.id}.xln.bedgraph.gz")  , emit: bedgraph
    tuple val(meta), path("${meta.id}.xpm.bedgraph.gz")  , emit: xpm
    path  "versions.yml"                                 , emit: versions

    script:
    """
    export LC_ALL=C

    samtools view -H $bam | grep "@SQ" | cut -f2,3 | sed 's/SN://g' | sed 's/LN://g' > genome.tsv
    bedtools bamtobed -i $bam > dedup.bed

    bedtools shift -m 1 -p -1 -i dedup.bed -g genome.tsv > shifted.bed
    bedtools genomecov -dz -strand + -5 -i shifted.bed -g genome.tsv \
        | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "+"}' \
        > pos.bed

    bedtools genomecov -dz -strand - -5 -i shifted.bed -g genome.tsv \
        | awk '{OFS="\t"}{print \$1, \$2, \$2+1, ".", \$3, "-"}' \
        > neg.bed

    cat pos.bed neg.bed | sort -k1,1 -k2,2n | gzip > ${meta.id}.xln.bed.gz
    zcat ${meta.id}.xln.bed.gz \
        | awk '{OFS = "\t"}{if (\$6 == "+") {print \$1, \$2, \$3, \$5} else {print \$1, \$2, \$3, -\$5}}' \
        | gzip \
        > ${meta.id}.xln.bedgraph.gz
    
    total=\$(samtools view -F 2304 -c ${bam})
    zcat ${meta.id}.xln.bedgraph.gz \
        | awk -v total=\${total} '{\$4 *= 1000000 / total; print}' \
        | gzip \
        > ${meta.id}.xpm.bedgraph.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
