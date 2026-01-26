rule rgi_bwt:
    input:
        reads = "resources/reads/{sample}.renamed.fastq.gz"
    output:
        report_dir = directory("results/{sample}/rgi_bwt")
    conda:
        "../envs/rgi.yaml"
    threads: 20
    shell:
        r"""
        set -euo pipefail
        mkdir -p {output.report_dir}

        echo "==> Running RGI BWT on ONT reads (output in {output.report_dir})."

        rgi bwt \
            --read_one {input.reads} \
            --output_file {output.report_dir}/{wildcards.sample} \
            --threads {threads} \
            --local \
            --clean
        """

