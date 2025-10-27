rule medaka:
    input:
        reads="results/{sample}/clean/{sample}_noh.fq.gz",
        ref="reference/{ref}.fasta"
    output:
        consensus="results/{sample}/medaka/{ref}/consensus.fasta"
    conda:
        "../envs/medaka.yaml"
    shell:
        "medaka_consensus -i {input.reads} -d {input.ref} -o results/{wildcards.sample}/medaka/{wildcards.ref} -m r1041_e82_400bps_fast_variant_g632 -t 8"
