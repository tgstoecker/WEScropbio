### Mapping with BWA-mem

rule map_reads_bwa:
    input:
        reads=["results/trimmed/{sample}.R1.fastq", "results/trimmed/{sample}.R2.fastq"],
        index=multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        "results/mapped/{sample}.bwa.sorted.bam",
    log:
        "logs/bwa_mem/{sample}.log",
    params:
        index="resources/genome.fasta",
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA'",
        sort="none",
        sort_order="coordinate",
        sort_extra="",
#    conda:
#        "../envs/bwa-mem.yaml",
    threads: 8
    # default is coordiante sort
#    shell:
#        "bwa mem "
#        "-t {threads} "
#        "-R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' " 
#        " {params.index} {input.reads} -o {output}"
    wrapper:
        "0.74.0/bio/bwa/mem"

### Mapping with Novoalign

rule map_reads_novoalign:
    input:
        reads=["results/trimmed/{sample}.R1.fastq", "results/trimmed/{sample}.R2.fastq"],
        index="resources/genome.novoalign.idx"
    output:
        "results/mapped/{sample}.novoalign.sam"
    conda:
        "../envs/novoalign.yaml"
    threads: 24
    shell:
        "novoalign -c {threads} -d {input.index} -f {input.reads} -i 203,20 -o SAM > {output}"

# sam to sorted bam 
rule novoalign_sam_sort_bam:
    input:
        "results/mapped/{sample}.novoalign.sam"
    output:
        temp("results/mapped/{sample}.novoalign.noreadgroup.sorted.bam"),
    threads: 24
    conda: "../envs/samtools.yaml"
    shell:
        "samtools sort -@ {threads} -O BAM {input} -o {output}"


rule novoalign_index_bam:
    input:
        "results/mapped/{sample}.novoalign.noreadgroup.sorted.bam"
    output:
        temp("results/mapped/{sample}.novoalign.noreadgroup.sorted.bam.bai")
    threads: 24
    shell:
        "samtools index -@ {threads} {input}"


rule novoalign_add_read_groups:
    input:
        "results/mapped/{sample}.novoalign.noreadgroup.sorted.bam"
    output:
        "results/mapped/{sample}.novoalign.sorted.bam"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard AddOrReplaceReadGroups "
        "I={input} "
        "O={output} "
        "RGID=1 "
        "RGLB=lib1 "
        "RGPL=illumina "
        "RGPU=unit1 "
        "RGSM={wildcards.sample}"




### Marking duplicates with Picard ###
### For both mappers!

rule mark_duplicates:
    input:
        "results/mapped/{sample}.{aligner}.sorted.bam",
    output:
        bam="results/dedup/{sample}.{aligner}.sorted.bam",
        metrics="results/qc/dedup/{sample}.{aligner}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.{aligner}.log",
    params:
        config["params"]["picard"]["MarkDuplicates"],
    wrapper:
        "0.59.2/bio/picard/markduplicates"


rule mark_duplicates_indexes:
    input:
        "results/dedup/{sample}.{aligner}.sorted.bam"
    output:
        "results/dedup/{sample}.{aligner}.sorted.bam.bai"
    threads: 24
    shell:
        "samtools index -@ {threads} {input}"


### Applying GATK Base Qulaity Recalibration ###


rule recalibrate_base_qualities:
    input:
        bam="results/dedup/{sample}.{aligner}.sorted.bam",
        bai="results/dedup/{sample}.{aligner}.sorted.bam.bai",
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=temp("results/recal/{sample}.{aligner}.grp"),
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"],
        java_opts="",
    log:
        "logs/gatk/baserecalibrator/{sample}.{aligner}.log",
    threads: 8
    wrapper:
        "0.62.0/bio/gatk/baserecalibratorspark"


#ruleorder: apply_bqsr > bam_index


rule apply_bqsr:
    input:
        bam="results/dedup/{sample}.{aligner}.sorted.bam",
        bai="results/dedup/{sample}.{aligner}.sorted.bam.bai",
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        recal_table="results/recal/{sample}.{aligner}.grp",
    output:
        bam=protected("results/recal/{sample}.{aligner}.sorted.bam"),
        bai="results/recal/{sample}.{aligner}.sorted.bai",
    log:
        "logs/gatk/gatk_applybqsr/{sample}.{aligner}.log",
    params:
        extra=config["params"]["gatk"]["applyBQSR"],  # optional
        java_opts="",  # optional
    wrapper:
        "0.62.0/bio/gatk/applybqsr"
