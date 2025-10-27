#AMR Pipeline for Metagenomic analysis
#Author : Thibaut Armel Chérif GNIMADI
#Affiliation : CERFIG 


#################### RUN PIPELINE ##########################
snakemake --use-conda -j 22 --rerun-incomplete

#################### RUN DAG ###############################
snakemake --report results/summary/workflow_report.html
snakemake --dag | dot -Tpng > results/workflow_dag.png
#################### END  ##################################
