nextflow.enable.dsl = 2
//params.dataset_name = "230731IMT.test"
params.samples = file("data/" + params.dataset_name + "/metadata/sample.sheet.csv")
params.phred_threshold = "30"
params.clip_5_prime = "9"
params.clip_3_prime = "6"
params.min_read_length = "40"
params.trim_galore_args = "--illumina --gzip"

params.bismark_args = "--non_directional --bam"
params.dedup_args = "--bam --single"

params.met_ext_buffer = "15G"
params.met_ext_args = "--bedGraph --CX --gazillion --gzip  --single-end"

params.cov2cyto_seqtype = "--nome-seq --gc"
params.bamsort = "-c"

params.genome_dir = file("data/genomes/ambystoma_mexicanum/")
params.gtf_path = file(params.genome_dir + "AmexG_v6_DD.gtf")
params.fa_path = file(params.genome_dir + "AmexG_v6_DD.fa")
params.index_path = file(params.genome_dir + "AmexG_v6_DD.idx")

process run_trimgalore {
    tag "${sample_id}"
    label "plow"
    publishDir "data/" + params.dataset_name + "/raw/", mode: 'copy', pattern: "*.fq.gz", saveAs: {filename -> "${sample_id}.trimmed.fastq.gz"}
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.txt", saveAs: {filename -> "${sample_id}.trimmed.txt"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.trim_galore.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.trim_galore.out"}
    if (workflow.profile == 'kcs' || workflow.profile == 'asc') {
        conda 'bioconda::trim-galore bioconda::cutadapt'
    } else {
        container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    }

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*.fq.gz"), emit: fastq
    tuple val(sample_id), path("*.txt"), emit: report
    stdout emit: verbiage

    script:
    """
    trim_galore ${params.trim_galore_args} \
        --quality ${params.phred_threshold} --length ${params.min_read_length} \
        --clip_r1 ${params.clip_5_prime} --three_prime_clip_R1 ${params.clip_3_prime} \
        ${fastq}
    """
}

process run_fastqc {
    tag "${sample_id}"
    label "psingle"
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.html", saveAs: {filename -> "${sample_id}.trimmed.html"}
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.zip", saveAs: {filename -> "${sample_id}.trimmed.zip"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.fastqc.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.fastqc.out"}
    //container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    container "${ workflow.containerEngine == 'singularity' ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
        'registry.hub.docker.com/biocontainers/fastqc:v0.11.9_cv8' }"

    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*.html"), emit: report

    script:
    """
    fastqc -f fastq -q ${fastq}
    """
}

process run_bismark {
    tag "${sample_id}"
    label "phigh"
    publishDir "data/" + params.dataset_name + "/raw/", mode: 'copy', pattern: "*.bam", saveAs: { filename -> "${sample_id}.bismark.bam"}
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.txt", saveAs: { filename -> "${sample_id}.bismark.txt"}
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.html", saveAs: {filename -> "${sample_id}.trimmed.html"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.bismark.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.bismark.out"}
    if (workflow.profile == 'kcs' || workflow.profile == 'asc') {
        conda 'bioconda::bowtie2=2.5.2  bioconda::samtools'
    } else {
        container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    }


    input:
    tuple val(sample_id), path(fastq)

    output:
    tuple val(sample_id), path("*.bam"), emit: bam
    tuple val(sample_id), path("*.txt"), emit: report

    script:
    """
    echo ${workflow.launchDir}
    n_threads=\$((${task.cpus}/2))
    bismark ${params.bismark_args} --temp_dir ./ \
        --output_dir ./ --genome ${params.genome_dir} \
        -p 8 --parallel 1 ${fastq}
    """
}

process run_dedup {
    tag "${sample_id}"
    label "psingle"
    publishDir "data/" + params.dataset_name + "/raw/", mode: 'copy', pattern: "*.bam", saveAs: { filename -> "${sample_id}.bismark.deduplicated.bam"}
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.txt", saveAs: { filename -> "${sample_id}.bismark.deduplicated.txt"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.deduplication.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.deduplication.out"}
    //container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    if (workflow.profile == 'kcs' || workflow.profile == 'asc') {
        conda 'bioconda::bowtie2  bioconda::samtools'
    } else {
        container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("*.bam"), emit: bam
    tuple val(sample_id), path("*.txt"), emit: report

    script:
    """
    deduplicate_bismark ${params.dedup_args} ${input_bam}
    """
}

process run_methlation_extract {
    tag "${sample_id}"
    label "pmedium"
    publishDir "data/" + params.dataset_name + "/raw/", mode: 'copy', pattern: "*.cov.gz", saveAs:  { filename -> "${sample_id}.bismark.deduplicated.cov.gz"}
    publishDir "data/" + params.dataset_name + "/raw/", mode: 'copy', pattern: "*.bedGraph.gz", saveAs: { filename -> "${sample_id}.bismark.deduplicated.bedGraph.gz"}
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*splitting_report.txt", saveAs: { filename -> "${sample_id}.bismark.deduplicated.splitting.txt"}
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*M-bias.txt", saveAs: { filename -> "${sample_id}.bismark.deduplicated.mbias.txt"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.met_extract.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.met_extract.out"}
    if (workflow.profile == 'kcs' || workflow.profile == 'asc') {
        conda 'bioconda::bowtie2=2.5.2  bioconda::samtools'
    } else {
        container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("*.cov.gz"), emit: cov
    tuple val(sample_id), path("*splitting_report.txt"), emit: report_splitting
    tuple val(sample_id), path("*M-bias.txt"), emit: report_mbias

    script:
    """
    bismark_methylation_extractor  --output ./ \
        --genome_folder ${params.genome_dir} --buffer ${params.met_ext_buffer} \
        --parallel ${task.cpus} ${params.met_ext_args} ${input_bam}
    """

}

process run_cov2cyto {
    tag "${sample_id}"
    label "pmedium"
    publishDir "data/" + params.dataset_name + "/cooked/", mode: 'copy', pattern: "*GpC.cov"
    publishDir "data/" + params.dataset_name + "/cooked/", mode: 'copy', pattern: "*CpG.cov"
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.txt"
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.cov2cyto.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.cov2cyto.out"}
    if (workflow.profile == 'kcs' || workflow.profile == 'asc') {
        conda 'bioconda::bowtie2=2.5.2  bioconda::samtools'
    } else {
        container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("*GpC.cov"), emit: cov_gpc
    tuple val(sample_id), path("*CpG.cov"), emit: cov_cpg
    tuple val(sample_id), path("*CpG_report.txt"), emit: report_cpg
    tuple val(sample_id), path("*GpC_report.txt"), emit: report_gpc

    script:
    """
    coverage2cytosine --genome_folder ${params.genome_dir} --dir ./  \
        -o ${sample_id} ${params.cov2cyto_seqtype} ${input_bam}
    """
}


process get_library_complexity {
    tag "${sample_id}"
    label "plow"
    publishDir "data/" + params.dataset_name + "/cooked/", mode: 'copy', pattern: "*.csv", saveAs: { filename -> "${sample_id}.complexity.csv"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.complexity.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.complexity.out"}
    if (workflow.profile == 'kcs' || workflow.profile == 'asc') {
        conda 'bioconda::bowtie2=2.5.2  bioconda::samtools bioconda::bedtools bioconda::preseq'
    } else {
        container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    }

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("*.csv"), emit: complexity

    script:
    """
    set +o pipefail;
    samtools sort ${input_bam} > ${sample_id}.bismark.sorted.bam
    bedtools bamtobed -i ${sample_id}.bismark.sorted.bam > ${sample_id}.bismark.sorted.bed
    preseq lc_extrap -o ${sample_id}.complexity.csv ${sample_id}.bismark.sorted.bed
    """

}

process run_multiqc {
    label "psingle"
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*.html"
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'copy', pattern: "*_data/"
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.multiqc.log"}
    publishDir "data/" + params.dataset_name + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.multiqc.out"}
    container 'quay.io/biocontainers/multiqc:1.18--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(input_report)
    val name

    output:
    path "*.html"

    script:
    """
    echo "input: ${input_dir}"
    multiqc --outdir ./ --filename ${params.dataset_name} --verbose --force \
        -c config/multiqc_config.yaml --interactive \
        --sample-names data/${params.dataset_name}/metadata/renaming.csv data/${params.dataset_name}/results/reports/
    """

}


workflow {
    Channel
        .fromPath( params.samples )
        .splitCsv( header: true, sep: ',' )
        .map { row -> tuple( row.sample, file(row.fastq) ) }
        .set { input }
    //input.view()

    run_trimgalore(input)
    fastq_trimmed = run_trimgalore.out.fastq
    trim_report = run_trimgalore.out.report
    //run_trimgalore.out.verbiage.view()

    run_fastqc(fastq_trimmed)
    fastq_report = run_fastqc.out.report

    run_bismark(fastq_trimmed)
    bismark_bam = run_bismark.out.bam
    bismark_report = run_bismark.out.report

    get_library_complexity(bismark_bam)

    run_dedup(bismark_bam)
    dedup_report = run_dedup.out.report
    dedup_bam = run_dedup.out.bam

    run_methlation_extract(dedup_bam)
    methylation_extract_cov = run_methlation_extract.out.cov
    methylation_extract_mbias = run_methlation_extract.out.report_mbias
    methylation_extract_splitting = run_methlation_extract.out.report_splitting

    run_cov2cyto(methylation_extract_cov)
    cov2cyto_report_cpg = run_cov2cyto.out.report_cpg
    cov2cyto_report_gpc = run_cov2cyto.out.report_gpc

    all_report = fastq_report.join(dedup_report)
        .join(trim_report)
        .join(bismark_report)
        .join(cov2cyto_report_cpg)
        .join(cov2cyto_report_gpc)
        .join(methylation_extract_mbias)
        .join(methylation_extract_splitting)
    run_multiqc(all_report)
}

