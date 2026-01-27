#AMR Pipeline for Metagenomic analysis
#Author : Thibaut Armel Chérif GNIMADI
#Affiliation : CERFIG 


#################### RUN PIPELINE ##########################
snakemake --until aggregate_rgi_samples  --cores 22 --use-conda

#################### RUN DAG ###############################
mkdir -p results/summary results
snakemake --snakefile workflow/Snakefile --report results/summary/workflow_report.html
snakemake --snakefile workflow/Snakefile --dag | dot -Tpng > results/workflow_dag.png
#################### END  ##################################
