rule samtools_merge_final_bams:
    input:
        expand(["results/recal/{sample}.{aligner}.sorted.bam"], sample=SAMPLES, aligner=ALIGNERS)
    output:
        "results/recal/merged/merged_recal.bam"
    params:
        "-c -p" # optional additional parameters as string
    threads:  # Samtools takes additional threads through its option -@
        24     # This value - 1 will be sent to -@
    wrapper:
        "0.74.0/bio/samtools/merge"


rule samtools_index:
    input:
        "results/recal/merged/merged_recal.bam"
    output:
        "results/recal/merged/merged_recal.bam.bai"
    params:
        "" # optional params string
    wrapper:
        "0.74.0/bio/samtools/index"


rule varlociraptor_preprocess:
    input:
        ref="resources/genome.fasta",
        ref_idx="resources/genome.fasta.fai",
        candidates="results/candidates_norm/all_norm.bcf",
        bam="results/recal/merged/merged_recal.bam",
        bai="results/recal/merged/merged_recal.bam.bai",
    output:
        "results/observations/observations.bcf",
    params:
#        extra=config["params"]["varlociraptor"]["preprocess"],
    log:
        "logs/varlociraptor/preprocess/preprocess.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --output {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs="results/observations/observations.bcf",
        scenario="config/scenario.yaml",
    output:
        "results/varlociraptor_calls/calls.bcf"
    log:
        "logs/varlociraptor/call/calls.log",
    params:
#        obs=lambda w, input: [
#            "{}={}".format(s, f) for s, f in zip(get_group_aliases(w), input.obs)
#        ],
        extra=config["params"]["varlociraptor_call"],
    conda:
        "../envs/varlociraptor.yaml"
#    benchmark:
#        "benchmarks/varlociraptor/call/{group}.{caller}.{scatteritem}.tsv"
    shell:
        "varlociraptor "
        "call variants generic --obs NA12878={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"

