#!/usr/bin/env nextflow

int threads = Runtime.getRuntime().availableProcessors()

Channel
    .fromPath(params.fasta)
    .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
    .into { fasta_to_index; fasta_minimap2; fasta_clairvoyante }
if(params.fai){
Channel
    .fromPath(params.fai)
    .ifEmpty{exit 1, "FASTA index file not found: ${params.fai}"}
    .set(fai)
}

// get relative path for model for clairvoyante
model = params.model.substring(params.model.lastIndexOf("/")+1)

// set model
model_data = "${params.model}.data-00000-of-00001"
model_index = "${params.model}.index"
model_meta = "${params.model}.meta"

Channel
    .fromPath(model_data)
    .ifEmpty { exit 1, "Model data file not found: ${model_data}" }
    .set { model_data }
Channel
    .fromPath(model_index)
    .ifEmpty { exit 1, "Model index file not found: ${model_index}" }
    .set { model_index }
Channel
    .fromPath(model_meta)
    .ifEmpty { exit 1, "Model meta file not found: ${model_meta}" }
    .set { model_meta }

Channel
      .fromPath(params.reads)
      .map { file -> tuple(file.baseName, file) }
      .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes" }
      .set { reads_minimap2 }

minimap2 = reads_minimap2.combine(fasta_minimap2)

if(!params.fai) {
  process preprocess_fai {
      tag "${fasta}"
      container 'lifebitai/samtools:latest'

      input:
      file(fasta) from fasta_to_index

      output:
      file("${fasta}.fai") into fai

      script:
      """
      samtools faidx $fasta
      """
  }
}

process minimap2 {
    tag "$reads"
    container 'evolbioinfo/minimap2:v2.14'

    cpus threads

    input:
    set val(name), file(reads), file(fasta) from minimap2

    output:
    set val(name), file("${name}.sam") into mapped_reads

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} $fasta $reads > ${name}.sam
    """
}

process bwa_sort {
    tag "$sam"
    container 'lifebitai/samtools:latest'

    input:
    set val(name), file(sam) from mapped_reads

    output:
    set val(name), file("${name}-sorted.bam") into sorted_bam

    """
    samtools sort -o ${name}-sorted.bam -O BAM $sam
    """
}

process mark_duplicates {
    tag "$bam"
    publishDir "${params.outdir}/marked_dup_bam", mode: 'copy'
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(bam) from sorted_bam

    output:
    set val(name), file("${name}.bam"), file("${name}.bai") into marked_bam_clairvoyante, marked_bam_sniffles
    file ("${name}.bam.metrics") into mark_dup_report

    """
    gatk MarkDuplicates \
    -I ${bam} \
    --CREATE_INDEX true \
    -M ${name}.bam.metrics \
    -O ${name}.bam
    """
}

clairvoyante = marked_bam_clairvoyante.merge(fasta_clairvoyante, fai, model_data, model_index, model_meta)

process clairvoyante {
    tag "$bam"
    publishDir "${params.outdir}/clairvoyante", mode: 'copy'
    container 'lifebitai/clairvoyante:latest'

    cpus threads

    input:
    set val(name), file(bam), file(bai), file(fasta), file(fai), file(model_data), file(model_index), file(model_meta) from clairvoyante

    output:
    set file("${name}.vcf.gz"), file("${name}.vcf.gz.tbi") into clairvoyante_vcf 

    // TODO: add optional param for `--bed_fn <file.bed> \`
    """
    clairvoyante.py callVarBamParallel \
       --chkpnt_fn $model \
       --ref_fn $fasta \
       --bam_fn $bam \
       --sampleName $name \
       --output_prefix $name \
       --threshold 0.125 \
       --minCoverage 4 \
       --tensorflowThreads ${task.cpus} \
       > commands.sh
    export CUDA_VISIBLE_DEVICES=""
    cat commands.sh | parallel -j${task.cpus}
    vcfcat ${name}*.vcf | vcfstreamsort | bgziptabix ${name}.vcf.gz
    """
}

process sniffles {
    tag "$bam"
    publishDir "${params.outdir}/sniffles", mode: 'copy'
    container 'lifebitai/sniffles:latest'

    cpus threads

    input:
    set val(name), file(bam), file(bai) from marked_bam_sniffles

    output:
    file("${name}.vcf") into sniffles_vcf 

    """
    sniffles --mapped_reads $bam --vcf ${name}.vcf --threads ${task.cpus}
    """
}