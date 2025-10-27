rule nanoplot:
    input:
        "resources/reads/{sample}.fastq.gz"
    output:
        directory("results/{sample}/qc/{sample}_Nanoplot")
    threads: 6
    conda:
        "../envs/nanoplot.yaml"
    log:
        "results/{sample}/qc/{sample}_Nanoplot/nanoplot.log"
    shell:
        r"""
        set -euo pipefail

        outdir="{output}"
        fq="{input}"
        logf="{log}"

        rm -rf "$outdir"
        mkdir -p "$outdir/tmp" "$outdir/.mplconfig"

        export TMPDIR="$outdir/tmp"
        export MPLBACKEND=Agg
        export MPLCONFIGDIR="$outdir/.mplconfig"

        # Exécution + redirection vers le log (pas d'option --logfile)
        set +e
        NanoPlot \
          --fastq "$fq" \
          -o "$outdir" \
          --threads {threads} \
          --huge \
          -f png \
          --verbose \
        > "$logf" 2>&1
        rc=$?
        set -e

        # Fallback: garder le pipeline vivant si NanoPlot échoue
        if [ $rc -ne 0 ]; then
            echo "[nanoplot] exit code $rc — generating placeholder" | tee -a "$logf"
            cat > "$outdir/index.html" <<HTML
<!doctype html><meta charset="utf-8">
<h3>NanoPlot — {wildcards.sample}</h3>
<p>Execution failed; keeping pipeline alive with a placeholder.</p>
<pre>
$(tail -n 200 "$logf" || true)
</pre>
HTML
        fi
        """
