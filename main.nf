#!/usr/bin/env nextflow

strandedness = params.forward_stranded ? "forward": "reverse"

log.info " Screen Preprocessing "
log.info "======================"
log.info "readDir                  : ${params.readDir}"
log.info "resultsDir               : ${params.resultsDir}"
log.info "library                  : ${params.library}"
log.info "demux barcode length     : ${params.barcode_demux_length}"
log.info "demux barcode mismatches : ${params.barcode_demux_mismatches}"
log.info "random barcode length    : ${params.barcode_random_length}"
log.info "spacer length            : ${params.spacer_length}"
log.info "strandedness             : ${strandedness}"
log.info "======================"

Channel
    .fromPath( "${params.readDir}/*.bam" )
    .map { file -> tuple( file.baseName, file ) }
    .set { bamFiles }

Channel
    .fromPath(params.library)
    .into { processLibraryFile ; mageckLibraryFile }

process process_library {

    tag "${library.baseName}"

    publishDir path: "${params.resultsDir}/library", mode: 'copy'

    // module params.R

    input:
    file(library) from processLibraryFile

    output:
    file("library.saf") into librarySafFile
    file("library.fasta") into libraryFastaFile

    script:
    """
    process_lib.R ${library} 'G'
    """
}

process build_bowtie_index {

    tag "${library_fasta.baseName}"

    // module params.bowtie2

    input:
    file(library_fasta) from libraryFastaFile

    output:
    file("bt2") into bt2Index

    script:
    """
    mkdir -p bt2
    bowtie2-build ${library_fasta} bt2/
    """
}

process bam_to_fastq {

    tag "${id}"

    // module params.samtools

    input:
    set val(id), file(bam) from bamFiles

    output:
    set val("${id}"), file("${id}.fastq") into fastqFiles

    script:
    """
    samtools fastq ${bam} > ${id}.fastq
    """
}

process trim_random_barcode {

    tag "${id}"

    // module params.fastx_toolkit

    input:
    set val(id), file(fastq) from fastqFiles

    output:
    set val("${id}"), file("${id}_random_barcode_trimmed.fastq") into randomBarcodeTrimmedFiles

    script:
    position = params.forward_stranded ? params.barcode_random_length + 1 : params.barcode_random_length
    flag_strandedness = params.forward_stranded ? "-f ${position}" : '-t {position}'

    """
    fastx_trimmer \
        ${flag_strandedness} \
        -Q33 \
        -i ${fastq} \
        -o ${id}_random_barcode_trimmed.fastq
    """
}

Channel
    .fromPath( "${params.readDir}/*.txt" )
    .map { file -> tuple(file.baseName, file) }
    .mix(randomBarcodeTrimmedFiles)
    .groupTuple()
    .set { demuxFiles }

process demultiplex_reads {

    tag "${id}"

    // module params.fastx_toolkit

    input:
    set val(id), file(files) from demuxFiles

    output:
    set val("${id}"), file('*.fastq') into splitFiles

    script:
    flag_strandedness = params.forward_stranded ? '--bol' : '--eol'
    num_mismatches = params.barcode_demux_mismatches

    """
    cat ${files[1]} | fastx_barcode_splitter.pl \
        --bcfile ${files[0]} \
        --prefix '' \
        --suffix .fastq \
        --mismatches ${num_mismatches} \
        ${flag_strandedness}
    """
}

def ungroupTuple = {
    def result = []
    def name = it[0]
    it[1].each { result << [name, it] }
    return result
 }

splitFiles
  .flatMap { it -> ungroupTuple(it) }
  .filter { it[1].baseName =~ /^(?!unmatched).*$/ }
  .set { flattenedSplitFiles }

process trim_barcode_and_spacer {

    tag "${id}_${fastq.baseName}"

    // module params.fastx_toolkit

    input:
    set val(id), file(fastq) from flattenedSplitFiles

    output:
    set val("${id}"), file("${id}_${fastq.baseName}.fastq") into spacerTrimmedFiles

    script:
    barcode_spacer_length = params.spacer_length + params.barcode_demux_length
    position = params.forward_stranded ? barcode_spacer_length + 1 : barcode_spacer_length
    flag_strandedness = params.forward_stranded ? "-f ${position}" : '-t {position}'

    """
    fastx_trimmer \
        ${flag_strandedness} \
        -Q33 \
        -i ${fastq} \
        -o ${id}_${fastq.baseName}.fastq
    """
}

spacerTrimmedFiles
    .into { alignmentSplitFiles ; compressSplitFiles }

process compress_reads {

    tag "${fastq.baseName}"

    publishDir path: "${params.resultsDir}/${id}/demultiplex", mode: 'copy'

    input:
    set val(id), file(fastq) from compressSplitFiles

    output:
    file("${fastq}.gz") into processedFiles

    script:
    """
    gzip -c ${fastq} > ${fastq}.gz
    """
}

process align_reads {

    tag "${fastq.baseName}"

    // module params.bowtie2
    // module params.samtools
    // module params.bamtools

    input:
    set val(id), file(fastq) from alignmentSplitFiles
    each file(index) from bt2Index

    output:
    set val("${id}"), file("${fastq.baseName}.bam") into alignedFiles

    script:
    """
    bowtie2 \
        -x ${index}/ \
        -L 18 \
        -N 0 \
        --seed 42 \
        ${fastq} \
    | samtools view -b - \
    | bamtools filter -tag AS:i:0 -out ${fastq.baseName}.bam

    """
}

alignedFiles
    .groupTuple()
    .set { groupedAlignedFiles }


process count_reads {

    tag "${id}"

    // module params.subread

    input:
    set val(id), file(bams) from groupedAlignedFiles
    each file(saf) from librarySafFile

    output:
    set val("${id}"), file("${id}_fc.txt") into countedFiles

    script:
    """
    featureCounts \
        -a ${saf} \
        -F SAF \
        -o ${id}_fc.txt \
        ${bams}
    """
}

process counts_to_mageck {

    tag "${id}"

    publishDir path: "${params.resultsDir}/${id}/counts", mode: 'copy', overwrite: 'true'

    // module params.R

    input:
    set val(id), file(counts) from countedFiles
    each file(library) from mageckLibraryFile

    output:
    file("${id}.txt") into mageckFiles

    script:
    """
    process_counts.R ${counts} ${library}
    """
}

process combine_counts {

    tag "all"

    publishDir path: "${params.resultsDir}/counts", mode: 'copy', overwrite: 'true'

    // module params.R

    input:
    file(counts) from mageckFiles.collect()

    output:
    file("counts.txt") into combinedCountsFile

    script:
    """
    combine_counts.R ${counts}
    """
}

workflow.onComplete {
	println ( workflow.success ? "DONE!" : "FAILED" )
}
