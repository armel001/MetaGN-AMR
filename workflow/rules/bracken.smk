rule bracken:
    input:
        kraken_report="results/{sample}/kraken/{sample}_report.txt"
    output:
        "results/{sample}/bracken/{sample}_bracken.txt"
    conda:
        "../envs/bracken.yaml"
    shell:
        "bracken -d databases/kraken2_8g_standard -i {input.kraken_report} -o {output} -r 100 -l S"