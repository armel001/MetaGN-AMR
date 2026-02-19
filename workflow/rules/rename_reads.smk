rule rename_reads:
    input:
        reads = "results/{sample}/clean/{sample}_noh.fq.gz"
    output:
        renamed = temp("results/{sample}/clean/{sample}.renamed.fastq.gz")
    log:
        "logs/rename_reads/{sample}.log"
    benchmark:
        "benchmarks/rename_reads/{sample}.txt"
    threads: 4
    resources:
        mem_mb  = 4000,    # 4 Gb RAM suffisant pour awk
        disk_mb = 50000    # 50 Gb espace temporaire
    conda:
        "../envs/python_filter.yaml"
    shell:
        """
        echo "[$(date)] Starting read renaming for {wildcards.sample}" > {log}
        echo "[$(date)] Input: {input.reads}" >> {log}

        # Compter les reads avant renommage
        n_reads=$(zcat {input.reads} | awk 'NR%4==1' | wc -l)
        echo "[$(date)] Reads to rename: $n_reads" >> {log}

        # Renommer les reads : @{sample}_1, @{sample}_2, ...
        # NR%4==1 : ligne header (@)
        # NR%4==2 : séquence
        # NR%4==3 : ligne + 
        # NR%4==4 : qualité
        zcat {input.reads} \
            | awk \
                'NR%4==1 {{printf("@{wildcards.sample}_%d\\n", ++i)}} \
                 NR%4!=1 {{print}}' \
            | bgzip -l 9 -@ {threads} \
            > {output.renamed} \
            2>> {log}

        # Vérifier que le nombre de reads est identique
        n_renamed=$(zcat {output.renamed} | awk 'NR%4==1' | wc -l)
        echo "[$(date)] Reads after renaming: $n_renamed" >> {log}

        if [ "$n_reads" != "$n_renamed" ]; then
            echo "[$(date)] ERROR: Read count mismatch! $n_reads vs $n_renamed" >> {log}
            exit 1
        fi

        echo "[$(date)] Read renaming completed for {wildcards.sample}" >> {log}
        """
