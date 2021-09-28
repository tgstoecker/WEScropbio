rule filter_odds:
    input:
        "results/varlociraptor_calls/calls.bcf",
    output:
        "results/varlociraptor_filter/calls.filtered_odds.bcf",
    params:
        events=["HET", "HOM"],
    log:
        "logs/filter-calls/posterior_odds/calls_odds_filter.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds very-strong < {input} > {output} 2> {log}"


VARTYPES=["SNV", "INS", "DEL"]

rule control_fdr:
    input:
        "results/varlociraptor_filter/calls.filtered_odds.bcf",
    output:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf",
    log:
        "logs/control-fdr/calls.control_fdr_{vartype}.log",
#    params:
#        query=get_fdr_control_params,
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr --local {input} --var {wildcards.vartype} "
        "--events HET HOM --fdr 0.0001 > {output} 2> {log}"


rule bcftools_index_final_calls:
    input:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf"
    output:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf.csi"
    params:
        extra=""  # optional parameters for bcftools index
    wrapper:
        "0.74.0/bio/bcftools/index"


rule bcftools_concat_final_calls:
    input:
        calls=expand(["results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf"], vartype=VARTYPES),
        index=expand(["results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf.csi"], vartype=VARTYPES),
    output:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.bcf"
    params:
        "-Ob --threads 24 -a -d none"  # optional parameters for bcftools concat (except -o)
    wrapper:
        "0.74.0/bio/bcftools/concat"


rule targets_filtered_bcftools_concat_final_calls:
    input:
        calls=expand(["results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf"], vartype=VARTYPES),
        index=expand(["results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf.csi"], vartype=VARTYPES),
    output:
        "results/varlociraptor_filter/targets_filt_calls.filtered_odds.fdr_controlled_all_vartypes.bcf"
    params:
        "-Ob --threads 24 -a -d none -R targets/targets.bed"  # optional parameters for bcftools concat (except -o)
    wrapper:
        "0.74.0/bio/bcftools/concat"


rule bcf_to_vcf_1:
    input:
        bcf="results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.bcf"
    output:
        vcf="results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.vcf"
    params:
        ""  # optional parameters for bcftools view (except -o)
    log:
        "logs/final1_bcf2vcf.log"
    wrapper:
        "0.74.0/bio/bcftools/view"


rule bcf_to_vcf_2:
    input:
        bcf="results/varlociraptor_filter/targets_filt_calls.filtered_odds.fdr_controlled_all_vartypes.bcf"
    output:
        vcf="results/varlociraptor_filter/targets_filt_calls.filtered_odds.fdr_controlled_all_vartypes.vcf"
    params:
        ""  # optional parameters for bcftools view (except -o)
    log:
        "logs/final2_bcf2vcf.log"
    wrapper:
        "0.74.0/bio/bcftools/view"



rule norm_final_vcf_1:
    input:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.vcf"
    output:
        "results/varlociraptor_filter/norm_final.vcf"
    params:
        "-m+both -f resources/genome.fasta -c s -d both --threads 24"  # optional parameters for bcftools norm (except -o)
    wrapper:
        "0.74.0/bio/bcftools/norm"


rule norm_final_vcf_2:
    input:
        "results/varlociraptor_filter/targets_filt_calls.filtered_odds.fdr_controlled_all_vartypes.vcf"
    output:
        "results/varlociraptor_filter/targets_filt_norm_final.vcf"
    params:
        "-m+both -f resources/genome.fasta -c s -d both --threads 24"  # optional parameters for bcftools norm (except -o)
    wrapper:
        "0.74.0/bio/bcftools/norm"
