### Mapping with BWA-mem

rule map_reads_bwa:
    input:
        reads=["results/trimmed/{sample}.R1.fastq.gz", "results/trimmed/{sample}.R2.fastq.gz"]
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

#rule map_reads_novoalign:
#novoalign -d "$Re1".nix -f Data/"$r1" Data/"$r2" -i 203,20 -o SAM

# sam to sorted bam 


### Marking duplicates with Picard ###
### For both mappers!

rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam",
    output:
        bam=temp("results/dedup/{sample}.sorted.bam"),
        metrics="results/qc/dedup/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        config["params"]["picard"]["MarkDuplicates"],
    wrapper:
        "0.59.2/bio/picard/markduplicates"



### Applying GATK Base Qulaity Recalibration ###


rule recalibrate_base_qualities:
    input:
#        bam=get_recalibrate_quality_input,
#        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        known="resources/variation.noiupac.vcf.gz",
        tbi="resources/variation.noiupac.vcf.gz.tbi",
    output:
        recal_table=temp("results/recal/{sample}.grp"),
    params:
        extra=config["params"]["gatk"]["BaseRecalibrator"],
        java_opts="",
    log:
        "logs/gatk/baserecalibrator/{sample}.log",
    threads: 8
    wrapper:
        "0.62.0/bio/gatk/baserecalibratorspark"


#ruleorder: apply_bqsr > bam_index


rule apply_bqsr:
    input:
#        bam=get_recalibrate_quality_input,
#        bai=lambda w: get_recalibrate_quality_input(w, bai=True),
        ref="resources/genome.fasta",
        ref_dict="resources/genome.dict",
        ref_fai="resources/genome.fasta.fai",
        recal_table="results/recal/{sample}.grp",
    output:
        bam=protected("results/recal/{sample}.sorted.bam"),
        bai="results/recal/{sample}.sorted.bai",
    log:
        "logs/gatk/gatk_applybqsr/{sample}.log",
    params:
        extra=config["params"]["gatk"]["applyBQSR"],  # optional
        java_opts="",  # optional
    wrapper:
        "0.62.0/bio/gatk/applybqsr"
