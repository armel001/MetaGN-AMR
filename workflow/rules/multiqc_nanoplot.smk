rule multiqc_nanoplot:
    input:
        nanoplot_dirs = expand("results/{sample}/qc/{sample}_Nanoplot", sample=config["samples_id"])
    output:
        data_dir = directory("results/summary/")
    conda:
        "../envs/multiqc.yaml"
    shell:
        """
        mkdir -p results/summary
        multiqc {input.nanoplot_dirs} -o {output.data_dir}
        """
