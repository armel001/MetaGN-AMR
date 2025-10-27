rule kraken2:
    input:
        "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        report="results/{sample}/kraken/{sample}_report.txt",
        output="results/{sample}/kraken/{sample}_output.txt"
    params:
        db="databases/kraken2_8g_standard"
    conda:
        "../envs/kraken2.yaml"
    shell:
        "kraken2 --db {params.db} --report {output.report} --output {output.output} {input}"
