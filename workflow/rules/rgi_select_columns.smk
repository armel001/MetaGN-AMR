rule select_rgi_columns:
    input:
        filtered = "results/{sample}/rgi_preprocessing/{sample}_filtered.txt"
    output:
        selected = "results/{sample}/rgi_preprocessing/{sample}_selected.txt"
    params:
        columns = config["rgi"]["columns_to_keep"]
    log:
        "logs/rgi_select_columns/{sample}.log"
    conda:
        "../envs/rgi_preprocessing.yaml"
    threads: 1
    script:
        "../scripts/select_rgi_columns.py"
