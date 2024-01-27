
nextflow.enable.dsl = 2
params.genome_dir = "data/genomes/ambystoma_mexicanum/"
params.gtf_path = file(params.genome_dir + "AmexG_v6_DD_split.gtf")
params.fa_path = file(params.genome_dir + "AmexG_v6_DD_split.fa")
params.index_path = file(params.genome_dir + "AmexG_v6_DD.idx")

params.fastqc_path = file("apps/FastQC/fastqc")

process trim_fastq {
    publishDir "data/" + params.dataset_name + "/raw/", mode: 'copy'
    cpus 1
    executor 'local'
    conda 'cutadapt'

    input:
    tuple val(sample_id), path(reads)

    output:
    file("*.trimmed.fastq.gz")

    script:
    r1_out = "${reads[0].simpleName}.trimmed.fastq.gz"
    r2_out = "${reads[1].simpleName}.trimmed.fastq.gz"
    """
    cutadapt --minimum-length=35 --pair-filter=any -a nextera_i5=CTGTCTCTTATA \
        -a polyA=AAAAAAAAAA -a polyT=TTTTTTTTTT -A nextera_i5=CTGTCTCTTATA \
        -A polyA=AAAAAAAAAA -A polyT=TTTTTTTTTT --times=2 \
        -o ${r1_out} -p ${r2_out} ${reads[0]} ${reads[1]}
    """
}

process run_fastqc {
    publishDir "data/" + params.dataset_name + "/results/reports/", mode: 'move'
    cpus 1
    executor 'local'

    input:
    tuple val(r1), val(r2)

    output:
    path '*.zip'
    path '*.html'
    val params.dataset_name

    script:
    """
    ${params.fastqc_path} ${r1} -q -t 1 --outdir ./
    ${params.fastqc_path} ${r2} -q -t 1 --outdir ./
    """
}


process align_to_genome_star {
    beforeScript 'export PATH=$PATH:bin/'
    publishDir "data/" + params.dataset_name + "/raw/", mode: 'copy'
    cpus 1
    memory '100 GB'

    input:
    tuple val(r1), val(r2)

    script:
    sample_id = ("${r1.simpleName}" - ~ /\w{3}\b/)
    """
    ulimit -v 900000000000
    STAR --genomeDir=${params.genome_dir} --readFilesCommand=zcat \
        --outSAMtype=SAM --outSAMattributes=RG NH HI AS NM nM MD MC jM jI \
        --outSAMunmapped=Within KeepPairs --runThreadN=16 --outSAMmultNmax=1 \
        --sjdbGTFfile=${params.gtf_path} --sjdbOverhang=100 \
        --outSAMattrRGline=ID:${sample_id}\
        --readFilesIn=${r1} ${r2}
    """
}

process align_kallisto {
    beforeScript 'export PATH=$PATH:bin/'
    publishDir "data/" + params.dataset_name + "/cooked/", mode: 'move'
    cpus 4
    memory '150 GB'
    executor 'slurm'

    input:
    tuple val(r1), val(r2)

    output:
    path '*'

    script:
    sample_id = ("${r1.simpleName}" - ~ /\w{3}\b/)
    """
    kallisto quant --bias -t 4 -b 100 \
        -i ${params.index_path} -t 4 -o ${sample_id} ${r1} ${r2}
    """
}

workflow {
    paired_raw_fastq = Channel.fromFilePairs("data/" + params.dataset_name + "/raw/L*_R{1,2}.fastq.gz", checkIfExists: true)
    trim_fastq(paired_raw_fastq) | (run_fastqc & align_kallisto)
}

workflow rerun_fastqc {
    trimmed_raw_fastq = Channel.fromFilePairs("data/" + params.dataset_name + "/raw/L*_R{1,2}.trimmed.fastq.gz", checkIfExists: true)
    run_fastqc(trimmed_raw_fastq)
}
