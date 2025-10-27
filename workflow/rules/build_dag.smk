rule workflow_summary:
    output:
        report="results/summary/workflow_report.html",
        dag="results/workflow_dag.png"
    conda:
        "../envs/dag.yaml"
    shell:
        """
        snakemake --snakefile workflow/Snakefile --report {output.report}
        snakemake --snakefile workflow/Snakefile --dag | dot -Tpng > {output.dag}
        """
