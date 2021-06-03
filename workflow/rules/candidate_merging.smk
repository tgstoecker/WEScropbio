SAMPLES = ["118711"]
ALIGNERS = ["bwa", "novoalign"]


rule concat_all_candidate_callsets:
    input:
        calls=expand(["results/candidate_calls/haplotypecaller/{sample}.{aligner}.bcf",
                      "results/candidate_calls/mpileup/{sample}.{aligner}.calls.bcf",
                      "results/candidate_calls/deep_variant/{sample}.{aligner}.reheader.bcf",
                      "results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.bcf",
                      "results/candidate_calls/delly/{sample}.{aligner}.delly.bcf"]
                     , sample=SAMPLES, aligner=ALIGNERS)
    output:
        "results/candidates_concat/all_concat.bcf"
    params:
        "-a --threads 24 -Ob"  # optional parameters for bcftools concat (except -o)
    wrapper:
        "0.74.0/bio/bcftools/concat"

