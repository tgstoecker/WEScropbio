rule concat_all_candidate_callsets:
    input:
        calls=expand(["results/candidate_calls/haplotypecaller/{sample}.{aligner}.reheader.bcf",
                      "results/candidate_calls/mpileup/{sample}.{aligner}.calls.reheader.bcf",
                      "results/candidate_calls/deep_variant/{sample}.{aligner}.reheader.bcf",
                      "results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.reheader.bcf",
                      "results/candidate_calls/delly/{sample}.{aligner}.delly.reheader.bcf"]
                     , sample=SAMPLES, aligner=ALIGNERS),
        index=expand(["results/candidate_calls/deep_variant/{sample}.{aligner}.reheader.bcf.csi",
                     "results/candidate_calls/haplotypecaller/{sample}.{aligner}.reheader.bcf.csi",
                     "results/candidate_calls/mpileup/{sample}.{aligner}.calls.reheader.bcf.csi",
                     "results/candidate_calls/freebayes/{sample}.{aligner}.freebayes.reheader.bcf.csi",
                     "results/candidate_calls/delly/{sample}.{aligner}.delly.reheader.bcf.csi"]
                     , sample=SAMPLES, aligner=ALIGNERS),
    output:
        "results/candidates_concat/all_concat.bcf"
    params:
        "-a --threads 24 -Ob"  # optional parameters for bcftools concat (except -o)
    wrapper:
        "0.74.0/bio/bcftools/concat"

