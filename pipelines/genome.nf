
params.genome_dir = file("data/genomes/ambystoma_mexicanum/")
params.gtf_path = file(params.genome_dir + "AmexG_v6_DD.gtf")
params.fa_path = file(params.genome_dir + "AmexG_v6_DD.fa")
params.tr_assembly_path = file(params.genome_dir + "AmexT_v47_dna.fa.zip")

process generate_genome_star {
    cpus = 12
    memory = "150 GB"

    script:
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 1000000000}" : ''

    """
    STAR --runMode genomeGenerate -—genomeDir ${params.genome_dir} --genomeFastaFiles ${params.fa_path} \
        --sjdbGTFfile ${params.gtf_path} --runThreadN ${task.cpus} ${memory}
    """
}

process generate_genome_kallisto {
    cpus = 12
    publishDir params.genome_dir, mode: 'move'
    memory '100 GB'
    executor 'slurm'

    output:
    file '*.idx'

    script:
    """
    kallisto index -i ${params.fa_path.simpleName}.idx ${params.tr_assembly_path}
    """
}

process genereate_genome_bismark {
    publishDir params.genome_dir + "/logs/", mode: 'copy', pattern: ".command.log", saveAs: {filename -> "${sample_id}.generate_genome.log"}
    publishDir params.genome_dir + "/logs/", mode: 'copy', pattern: ".command.out", saveAs: {filename -> "${sample_id}.generate_genome.out"}
    //container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    if (workflow.profile == 'kcs' || workflow.profile == 'asc') {
        conda 'bioconda::bowtie2  bioconda::samtools'
    } else {
        container 'docker.io/onurcanbektas/axolotl.rulands.yun:latest'
    }

    publishDir params.genome_dir, mode: 'move'
    memory '200GB'
    cpus '32'
    time { 48.hour * task.attempt }
    executor 'slurm'

    output:
    file '.command.log'
    file '.command.out'

    script:
    """
    bismark_genome_preparation --bowtie2 --verbose ${params.genome_dir} --large-index `expr ${task.cpus} / 2`
    """
}

workflow {
    //generate_genome_kallisto()
    genereate_genome_bismark()
}
