# Example:
# cd /path/to/project_A/analysis/dataset1/

# Run snakemake, specifying the location of the Snakefile and using the local config.yaml
snakemake --snakefile /home/lag36/lag_paulsson-lab/nanopore-env_code/library_nanopore_analysis/Snakefile \
          --configfile config.yaml \
          --profile /home/lag36/lag_paulsson-lab/nanopore-env_code/library_nanopore_analysis/slurm/ 
#          -n
