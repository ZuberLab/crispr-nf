#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ================================================================
     crispr-nf
    ================================================================

    Pre-processing of CRISPR-Cas9 / shRNA screening data

    Usage:
    nextflow run ZuberLab/crispr-nf

    Options:
        --readDir               Directory containing the raw reads. Valid input
                                files are BAM and compressed FASTQ. (Defaults to
                                'reads')

        --library               Path to sgRNA / shRNA library file. (Defaults to
                                'library.txt')
                                The following columns are required:
                                    - id: unique name of sgRNA / shRNA
                                    - gene: gene targeted by sgRNA / shRNA
                                    - sequence: nucleotide sequence of sgRNA / shRNA

        --barcodes              Path to file containing barcodes for demultiplexing.
                                (Defaults to 'barcodes.txt')
                                The following columns are required:
                                    - lane: name of BAM / FASTQ input file
                                    - sample_name: name of demultiplexed sample
                                    - barcode: nucleotide sequence of sample barcode

        --forward_stranded         Orientation of reads in BAM / FASTQ files.
                                   (Defaults to true)

        --barcode_random_length    Number of nucleotides in random barcode
                                   (Defaults to 6)

        --barcode_demux_mismatches Number of mismatches allowed during demultiplexing
                                   of barcode. (Defaults to 1)

        --barcode_demux_length     Number of nucleotides in sample barcode.
                                   (Defaults to 4)

        --spacer_length            Number of nucleotides in spacer sequence between
                                   barcodes and sgRNA / shRNA sequence. (Defaults to 20)

        --padding_base             Nucleotide used for padding if sgRNA / shRNA are of
                                   unequal length. Must be one of G, C, T, and A.
                                   (Defaults to G)

        --resultsDir               Directory to save results to. (Defaults to
                                   'results')

    Profiles:
        standard        local execution with singularity
        sge             SGE execution with singularity
        slurm           SLURM execution with singularity
        local           local execution without singularity

    Docker:
    zuberlab/crispr-nf:latest

    Author:
    Jesse J. Lipp (jesse.lipp@imp.ac.at)

    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}


strandedness = params.forward_stranded ? "forward": "reverse"

log.info " Screen Preprocessing "
log.info "======================"
log.info "readDir                  : ${params.readDir}"
log.info "resultsDir               : ${params.resultsDir}"
log.info "library file             : ${params.library}"
log.info "barcode file             : ${params.barcodes}"
log.info "demux barcode length     : ${params.barcode_demux_length}"
log.info "demux barcode mismatches : ${params.barcode_demux_mismatches}"
log.info "random barcode length    : ${params.barcode_random_length}"
log.info "spacer length            : ${params.spacer_length}"
log.info "strandedness             : ${strandedness}"
log.info "padding base             : ${params.padding_base}"
log.info "======================"

Channel
    .fromPath( "${params.readDir}/*.bam" )
    .map { file -> tuple( file.baseName, file ) }
    .set { bamInputFiles }

Channel
    .fromPath( "${params.readDir}/*.{fastq,fq}.gz" )
    .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
    .set { fastqInputFiles }

Channel
    .fromPath(params.barcodes)
    .into { processBarcodeFiles; combineBarcodeFiles }

Channel
    .fromPath(params.library)
    .into { processLibraryFile ; combineLibraryFile }

process process_library {

    tag "${library.baseName}"

    input:
    file(library) from processLibraryFile

    output:
    file("library.saf") into librarySafFile
    file("library.fasta") into libraryFastaFile

    script:
    padding_base = params.padding_base
    """
    process_library.R ${library} ${strandedness} ${padding_base}
    """
}

process process_barcodes {

    tag "${barcodes.baseName}"

    input:
    file(barcodes) from processBarcodeFiles

    output:
    file("*.txt") into demuxBarcodeFiles

    script:
    """
    process_barcodes.R ${barcodes}
    """
}

process bowtie_index {

    tag "${library_fasta.baseName}"

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

    tag { lane }

    input:
    set val(lane), file(bam) from bamInputFiles

    output:
    set val(lane), file("${lane}.fastq.gz") into fastqFilesFromBam

    script:
    """
    samtools fastq ${bam} | gzip -c > ${lane}.fastq.gz
    """
}

fastqFilesFromBam
    .mix(fastqInputFiles)
    .set { fastqFiles }

process trim_random_barcode {

    tag { lane }

    input:
    set val(lane), file(fastq) from fastqFiles

    output:
    set val(lane), file("${lane}_trimmed.fastq.gz") into randomBarcodeTrimmedFiles

    script:
    position = params.forward_stranded ? params.barcode_random_length + 1 : params.barcode_random_length
    flag_strandedness = params.forward_stranded ? "-f ${position}" : "-t ${position}"
    flag_strandedness = params.barcode_random_length > 0 ? flag_strandedness : '-f 1'
    
    """
    zcat ${fastq} | fastx_trimmer \
        ${flag_strandedness} \
        -Q33 \
        -z \
        -o ${lane}_trimmed.fastq.gz
    """
}

demuxBarcodeFiles
    .flatten()
    .map { file -> [file.baseName, file] }
    .concat(randomBarcodeTrimmedFiles)
    .groupTuple()
    .set { demuxFiles }

process demultiplex {

    tag "${lane}"

    input:
    set val(lane), file(files) from demuxFiles

    output:
    set val(lane), file('*.fq.gz') into splitFiles

    script:
    flag_strandedness = params.forward_stranded ? '--bol' : '--eol'
    num_mismatches = params.barcode_demux_mismatches

    """
    zcat ${files[1]} | fastx_barcode_splitter.pl \
        --bcfile ${files[0]} \
        --prefix ${lane}_ \
        --suffix .fq \
        --mismatches ${num_mismatches} \
        ${flag_strandedness}
    gzip *.fq
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
  .filter { it[1].baseName =~ /^(?!.*_unmatched).*$/ }
  .map { lane, file -> tuple(lane, file.name.replaceAll(/\.fq\.gz/, ''), file) }
  .into { flattenedSplitFiles; fastqcSplitFiles }


process trim_barcode_and_spacer {

    tag { id }

    input:
    set val(lane), val(id), file(fastq) from flattenedSplitFiles

    output:
    set val(lane), val(id), file("${id}.fastq.gz") into spacerTrimmedFiles

    script:
    barcode_spacer_length = params.spacer_length + params.barcode_demux_length
    position = params.forward_stranded ? barcode_spacer_length + 1 : barcode_spacer_length
    flag_strandedness = params.forward_stranded ? "-f ${position}" : "-t ${position}"

    """
    zcat ${fastq} | fastx_trimmer \
        ${flag_strandedness} \
        -Q33 \
        -z \
        -o ${id}.fastq.gz
    """
}

process align {

    tag { id }

    input:
    set val(lane), val(id), file(fastq) from spacerTrimmedFiles
    each file(index) from bt2Index

    output:
    set val(lane), val(id), file("${id}.bam") into alignedFiles
    file "${id}.log" into alignResults

    script:
    """
    bowtie2 \
        --threads ${task.cpus} \
        -x ${index}/ \
        -L 18 \
        -N 0 \
        --seed 42 \
        <(zcat ${fastq}) 2> ${id}.log \
    | samtools view -b - \
    | bamtools filter -tag AS:i:0 -out ${id}.bam

    """
}

alignedFiles
    .map { lane, id, file -> tuple(lane, file) }
    .groupTuple()
    .set { groupedAlignedFiles }

process count {

    tag { lane }

    publishDir path: "${params.resultsDir}/${lane}/counts",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(lane), file(bams) from groupedAlignedFiles
    each file(saf) from librarySafFile

    output:
    file("${lane}.txt") into countedFiles
    file("${lane}.txt.summary") into featureCountsResults

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -a ${saf} \
        -F SAF \
        -o ${lane}.txt \
        ${bams}
    """
}

process combine_counts {

    tag 'all'

    publishDir path: "${params.resultsDir}/counts", mode: 'copy', overwrite: 'true'

    input:
    file(counts) from countedFiles.collect()
    file(library) from combineLibraryFile
    file(barcodes) from combineBarcodeFiles

    output:
    file("counts.rds") into combinedCountsFile
    file("counts_mageck.txt") into combinedMageckFile

    script:
    """
    combine_counts.R ${library} ${barcodes} ${counts}
    """
}

fastqcSplitFiles
    .map { lane, id, file -> tuple(lane, file) }
    .groupTuple()
    .set { fastqcSplitFiles }

process fastqc {

    tag { lane }

    input:
    set val(lane), file(fastq) from fastqcSplitFiles

    output:
    file "*_fastqc.{zip,html}" into fastqcResults

    script:
    """
    fastqc -q ${fastq}
    """
}

process multiqc {

    tag 'all'

    publishDir "${params.resultsDir}/multiqc", mode: 'copy'

    input:
    file (fastqc: 'fastqc/*') from fastqcResults.collect()
    file (align: 'align/*') from alignResults.collect()
    file (featurecounts: 'featureCounts/*') from featureCountsResults.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"
    file '.command.err' into multiqc_stderr

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc -f -x *.run .
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
