rule norm_merged_candidate_callset:
    input:
        "results/candidates_concat/all_concat.bcf"
    output:
        "results/candidates_norm/all_norm.bcf"
    params:
        "--threads 24 -f resources/genome.fasta -d none -Ob -c s"  # optional parameters for bcftools norm (except -o)
    wrapper:
        "0.74.0/bio/bcftools/norm"
