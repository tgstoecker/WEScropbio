rule cutadapt_pe:
    input:
        ["rawreads/{sample}_1.fq.gz", "rawreads/{sample}_2.fq.gz"]
    output:
        fastq1="results/trimmed/{sample}.R1.fastq.gz",
        fastq2="results/trimmed/{sample}.R2.fastq.gz",
        qc="results/trimmed/{sample}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}.log",
    params:
        adapters="-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        extra="--minimum-length 1 -q 20"
    threads: 8
    wrapper:
        "0.74.0/bio/cutadapt/pe"

