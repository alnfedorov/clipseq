include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_RNALIB } from '../../modules/local/star_genomegenerate'
include { STAR_GENOMEGENERATE as STAR_GENOMEGENERATE_GENOME } from '../../modules/local/star_genomegenerate'

workflow PREPARE_GENOME {
    take:
    rnalib_fasta         // value: path(rnalib.fasta)
    rnalib_index         // value: path(rnalib STAR index)

    genome_fasta        // value: path(genome.fasta)
    genome_gtf          // value: path(genome.gtf)
    genome_index        // value: path(genome STAR index)
    
    main:
    versions = Channel.empty()

    if (!rnalib_index) {
        STAR_GENOMEGENERATE_RNALIB(
            file(rnalib_fasta, checkIfExists: true), 
            file('NA')
        )
        rnalib_index = STAR_GENOMEGENERATE_RNALIB.out.index
        versions = versions.mix(STAR_GENOMEGENERATE_RNALIB.out.versions)
    } else {
        rnalib_index = file(rnalib_index, checkIfExists: true)
    }

    if (!genome_index) {
        STAR_GENOMEGENERATE_GENOME(
            file(genome_fasta, checkIfExists: true), 
            file(genome_gtf, checkIfExists: true)
        )
        
        genome_index = STAR_GENOMEGENERATE_GENOME.out.index
        versions = versions.mix(STAR_GENOMEGENERATE_GENOME.out.versions)
    } else {
        genome_index = file(genome_index, checkIfExists: true)
    }

    emit:
    rnalib_index  = rnalib_index    // channel [ path(rnalib STAR index) ]
    genome_index = genome_index   // channel [ path(genome STAR index) ]

    versions = versions.ifEmpty(null) // channel: [ path(versions.yml) ]
}
