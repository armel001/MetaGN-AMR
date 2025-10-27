rule merge_rgi_results:
    input:
        expand("results/{sample}/rgi/{sample}.txt", sample=SAMPLES)
    output:
        merged = "results/summary/all_rgi_merged.tsv"
    conda:
        "../envs/python_plot.yaml"
    shell:
        r"""
        set -euo pipefail
        
        mkdir -p $(dirname {output.merged})

        echo "==> Extracting header from {input[0]}..."
        head -n 1 {input[0]} > {output.merged}

        echo "==> Merging content from all RGI files..."
        tail -n +2 {input} >> {output.merged}
        
        echo "Merged RGI table created: {output.merged}"
        """
