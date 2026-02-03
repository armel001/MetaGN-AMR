#AMR Pipeline for Metagenomic analysis
#Author : Thibaut Armel Chérif GNIMADI
#Affiliation : CERFIG 


#################### RUN PIPELINE ##########################
snakemake  --cores 22 --use-conda --rerun-triggers mtime

#################### RUN DAG ###############################
mkdir -p results/summary results
snakemake --snakefile workflow/Snakefile --report results/summary/workflow_report.html
snakemake --snakefile workflow/Snakefile --dag | dot -Tpng > results/workflow_dag.png
#################### END  ##################################
