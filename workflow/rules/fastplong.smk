rule fastplong:
    input:
        ancient("/data/armel/mythesis-data/data/{sample}.fastq.gz")
    output:
        fq   = temp("results/{sample}/fastplong/fp_trimmed.fastq.gz"),
        html = "results/{sample}/fastplong/fastplong.html",
        json = "results/{sample}/fastplong/fastplong.json"
    conda:
        "../envs/fastplong.yaml"
    threads: 16
    resources:
        mem_mb  = 16000,   # 16 Gb RAM
        disk_mb = 60000    # 60 Gb espace temporaire
    log:
        "logs/fastplong/{sample}.log"
    benchmark:
        "benchmarks/fastplong/{sample}.txt"
    shell:
        """
        echo "[$(date)] Starting fastplong for {wildcards.sample}" > {log}
        echo "[$(date)] Input: {input}" >> {log}

        fastplong \
            -i {input} \
            -o {output.fq} \
            -q 10 \
            -t {threads} \
            -h {output.html} \
            -j {output.json} \
            2>> {log}

        echo "[$(date)] Reads before filtering:" >> {log}
        echo "  $(zcat {input} | awk 'NR%4==1' | wc -l) reads" >> {log}

        echo "[$(date)] Reads after filtering:" >> {log}
        echo "  $(zcat {output.fq} | awk 'NR%4==1' | wc -l) reads" >> {log}

        echo "[$(date)] fastplong completed for {wildcards.sample}" >> {log}
        """
