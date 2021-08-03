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
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"


VARTYPES=["SNV", "INS", "DEL", "MNV", "BND", "INV", "DUP", "REP"]

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
        "--events HET HOM --fdr 0.05 > {output} 2> {log}"


rule bcftools_index_final_calls:
    input:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf"
    output:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf.csi"
    params:
        extra=""  # optional parameters for bcftools index
    wrapper:
        "0.74.0/bio/bcftools/index"


#Agilent target file
rule bcftools_concat_final_calls:
    input:
        calls=expand(["results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf"], vartype=VARTYPES),
        index=expand(["results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_{vartype}.bcf.csi"], vartype=VARTYPES),
    output:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.bcf"
    params:
        "-Ob --threads 24 -a -d none -R targets/targets.bed"  # optional parameters for bcftools concat (except -o)
    wrapper:
        "0.74.0/bio/bcftools/concat"


rule bcf_to_vcf:
    input:
        bcf="results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.bcf"
    output:
        vcf="results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.vcf"
    params:
        ""  # optional parameters for bcftools view (except -o)
    log:
        "logs/final_bcf2vcf.log"
    wrapper:
        "0.74.0/bio/bcftools/view"


rule norm_final_vcf:
    input:
        "results/varlociraptor_filter/calls.filtered_odds.fdr_controlled_all_vartypes.vcf"
    output:
        "results/varlociraptor_filter/norm_final.vcf"
    params:
#        "-c w --threads 24"  # optional parameters for bcftools norm (except -o)
        "-m+both -f resources/genome.fasta -c s -d both --threads 24"  # optional parameters for bcftools norm (except -o)
    wrapper:
        "0.74.0/bio/bcftools/norm"
