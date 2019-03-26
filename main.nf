#!/usr/bin/env nextflow

int threads = Runtime.getRuntime().availableProcessors()

Channel
    .fromPath(params.fasta)
    .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
    .into {fasta_to_index; fasta_minimap2; fasta_clairvoyante; fasta_sniffles}
if(params.fai){
Channel
    .fromPath(params.fai)
    .ifEmpty{exit 1, "FASTA index file not found: ${params.fai}"}
    .into {fai_clairvoyante; fai_sniffles}
}

// set the trained model files for Clairvoyante
params.data = params.model ? params.models[ params.model ].data ?: false : false
if (params.data) {
    Channel.fromPath(params.data)
           .ifEmpty { exit 1, "Trained model data file for Clairvoyante not found: ${params.data}" }
           .set { model_data }
}
params.index = params.model ? params.models[ params.model ].index ?: false : false
if (params.index) {
    Channel.fromPath(params.index)
           .ifEmpty { exit 1, "Trained model index file for Clairvoyante not found: ${params.index}" }
           .set { model_index }
}
params.meta = params.model ? params.models[ params.model ].meta ?: false : false
if (params.meta) {
    Channel.fromPath(params.meta)
           .ifEmpty { exit 1, "Trained model meta file for Clairvoyante not found: ${params.meta}" }
           .set { model_meta }
}

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
      file("${fasta}.fai") into (fai_clairvoyante, fai_sniffles)

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
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(bam) from sorted_bam

    output:
    set val(name), file("${name}-marked_dup.bam"), file("${name}-marked_dup.bai") into marked_bam_clairvoyante, marked_bam_sniffles, marked_bam_svim
    file "${name}.bam.metrics" into mark_dup_report

    """
    gatk MarkDuplicates \
    -I ${bam} \
    --CREATE_INDEX true \
    -M ${name}.bam.metrics \
    -O ${name}-marked_dup.bam
    """
}

clairvoyante = marked_bam_clairvoyante.merge(fasta_clairvoyante, fai_clairvoyante, model_data, model_index, model_meta)

process clairvoyante {
    tag "$bam"
    publishDir "${params.outdir}/clairvoyante", mode: 'copy'
    container 'lifebitai/clairvoyante:latest'

    cpus threads

    input:
    set val(name), file(bam), file(bai), file(fasta), file(fai), file(model_data), file(model_index), file(model_meta) from clairvoyante

    output:
    set file("${name}.vcf.gz"), file("${name}.vcf.gz.tbi") into clairvoyante_vcf 

    script:
    model = model_index.getName().substring(0, model_index.getName().lastIndexOf(".index"))
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

sniffles_preprocessing = marked_bam_sniffles.merge(fasta_sniffles, fai_sniffles)

// recompute the MD string which is required for sniffles
process sniffles_preprocessing {
    tag "$bam"
    publishDir "${params.outdir}/bam", mode: 'copy'
    container 'lifebitai/samtools:latest'

    input:
    set val(name), file(bam), file(bai), file(fasta), file(fai) from sniffles_preprocessing

    output:
    set val(name), file("${name}.bam"), file("${name}.bam.bai") into bam_md_sniffles 

    """
    samtools calmd $bam $fasta | samtools view -S -b - > ${name}.bam
    samtools index ${name}.bam
    """
}

process sniffles {
    tag "$bam"
    publishDir "${params.outdir}/sniffles", mode: 'copy'
    container 'lifebitai/sniffles:latest'

    cpus threads

    input:
    set val(name), file(bam), file(bai) from bam_md_sniffles

    output:
    file("sniffles_${name}.vcf") into sniffles_vcf

    """
    sniffles --mapped_reads $bam --vcf sniffles_${name}.vcf -s ${params.min_support} --threads ${task.cpus}
    """
}

sniffles_vcf
    .map { file -> tuple(file.baseName, file) }
    .into { sniffles_vcf_length; sniffles_vcf_carriers }

process svim {
    tag "$bam"
    container 'lifebitai/svim:latest'

    input:
    set val(name), file(bam), file(bai) from marked_bam_svim

    output:
    set val(name), file("${name}/final_results.vcf") into svim_vcf 

    """
    svim alignment ${name} $bam
    """
}

process filter_svim {
    tag "$vcf"
    publishDir "${params.outdir}/svim", mode: 'copy'
    container 'lifebitai/svim:latest'

    input:
    set val(name), file(vcf) from svim_vcf

    output:
    file("svim_${name}.vcf") into svim_filtered_vcf

    """
    cat $vcf | \
         awk '{{ if(\$1 ~ /^#/) {{ print \$0 }} \
         else {{ if(\$6>10) {{ print \$0 }} }} }}' > svim_${name}.vcf
    """
}

svim_filtered_vcf
    .map { file -> tuple(file.baseName, file) }
    .into { svim_filtered_vcf_length; svim_filtered_vcf_carriers}

process sv_length_plot {
    publishDir "${params.outdir}/plots", mode: 'copy'
    container 'lifebitai/sv-plots:latest'

    input:
    set val(sniffles_name), file(sniffles_vcf) from sniffles_vcf_length
    set val(svim_name), file(svim_vcf) from svim_filtered_vcf_length

    output:
    set file("SV-length_${sniffles_name}.png"), file("${sniffles_name}.txt"), file("SV-length_${svim_name}.png"), file("${svim_name}.txt") into sv_length_plots

    """
    SV-length-plot.py $sniffles_vcf -o SV-length_${sniffles_name}.png -c ${sniffles_name}.txt
    SV-length-plot.py $svim_vcf -o SV-length_${svim_name}.png -c ${svim_name}.txt
    """
}


process sv_carriers_plot {
    publishDir "${params.outdir}/plots", mode: 'copy'
    container 'lifebitai/sv-plots:latest'

    input:
    set val(sniffles_name), file(sniffles_vcf) from sniffles_vcf_carriers
    set val(svim_name), file(svim_vcf) from svim_filtered_vcf_carriers

    output:
    set file("SV-carriers_${sniffles_name}.png"), file("SV-carriers_${svim_name}.png") into sv_carriers_plots

    """
    SV-carriers-plot.py $sniffles_vcf -o SV-carriers_${sniffles_name}.png
    SV-carriers-plot.py $svim_vcf -o SV-carriers_${svim_name}.png
    """
}

process multiqc {
    tag "multiqc_report.html"

    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    container 'ewels/multiqc:v1.7'

    input:
    file bam_metrics from mark_dup_report

    when: 
    !params.skip_multiqc

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc . -m qualimap -m picard -m gatk -m bcftools
    """
}


process deploit_report {
    publishDir "${params.outdir}/Visualisations", mode: 'copy'

    container 'lifebitai/vizjson:latest'

    input:
    set file(sv_sniffles_length_plot), file(sniffles_table), file(sv_length_svim_plot), file(svim_table) from sv_length_plots
    set file(sniffles_carriers_plot), file(svim_carriers_plot) from sv_carriers_plots

    output:
    file '.report.json' into results

    when:
    !params.skip_deploit_report

    script:
    """
    sniffles_title=\$(head -n 1 $sniffles_table)
    svim_title=\$(head -n 1 $svim_table)

    tail -n +2 $sniffles_table > sniffles.tsv
    tail -n +2 $svim_table > svim.tsv

    tsv2csv.py < sniffles.tsv > sniffles.csv
    tsv2csv.py < svim.tsv > svim.csv

    csv2json.py sniffles.csv "\${sniffles_title}" "${sniffles_table}.json" False
    csv2json.py svim.csv "\${svim_title}" "${svim_table}.json" False

    img2json.py "results/plots/${sv_sniffles_length_plot}" "Sniffles Length Plot" "${sv_sniffles_length_plot}.json"
    img2json.py "results/plots/${sv_length_svim_plot}" "SVIM Length Plot" "${sv_length_svim_plot}.json"
    img2json.py "results/plots/${sniffles_carriers_plot}" "Sniffles Carriers Plot" "${sniffles_carriers_plot}.json"
    img2json.py "results/plots/${svim_carriers_plot}" "SVIM Carriers Plot" "${svim_carriers_plot}.json"

    combine_reports.py .
    """
}