rule gatk_calling:
    input:
        # single or list of bam files
        bam="mapped/{sample}.bam",
        ref="genome.fasta"
        # known="dbsnp.vcf"  # optional
    output:
        vcf="calls/{sample}.vcf",
    log:
        "logs/gatk/haplotypecaller/{sample}.log"
    params:
        extra="",  # optional
        java_opts="", # optional
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    wrapper:
        "0.74.0/bio/gatk/haplotypecaller"

rule samtools_mpileup:
    input:
        # single or list of bam files
        bam="mapped/{sample}.bam",
        reference_genome="genome.fasta"
    output:
        "mpileup/{sample}.mpileup.gz"
    log:
        "logs/samtools/mpileup/{sample}.log"
    params:
        extra="-d 10000",  # optional
    wrapper:
        "0.74.0/bio/samtools/mpileup"


rule deep_variant_calling:
    input:
        bam="mapped/{sample}.bam",
        ref="genome/genome.fasta"
    output:
        vcf="calls/{sample}.vcf.gz"
    params:
        model="wgs",   # {wgs, wes, pacbio, hybrid}
        sample_name=lambda w: w.sample, # optional
        extra=""
    threads: 2
    log:
        "logs/deepvariant/{sample}/stdout.log"
    wrapper:
        "0.74.0/bio/deepvariant"


rule freebayes:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        regions=get_regions(),
        # you can have a list of samples here
        samples=lambda w: get_group_bams(w),
        index=lambda w: get_group_bams(w, bai=True),
    output:
        "results/candidate-calls/{group}.freebayes.bcf",
    log:
        "logs/freebayes/{group}.log",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra="--pooled-continuous --min-alternate-count 1 --min-alternate-fraction {}".format(
            config["params"]["freebayes"].get("min_alternate_fraction", "0.05")
        ),
    threads: workflow.cores - 1  # use all available cores -1 (because of the pipe) for calling
    wrapper:
        "0.68.0/bio/freebayes"


rule delly:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        samples=lambda w: get_group_bams(w),
        index=lambda w: get_group_bams(w, bai=True),
        exclude=get_excluded_regions(),
    output:
        "results/candidate-calls/{group}.delly.bcf",
    log:
        "logs/delly/{group}.log",
    params:
        extra=config["params"].get("delly", ""),
    threads: lambda _, input: len(input.samples)  # delly parallelizes over the number of samples
    wrapper:
        "0.68.0/bio/delly"
