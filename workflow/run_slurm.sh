# SLURM specifications made in default.cluster.yaml & the individual rules
# GRB_LICENSE_FILE=/share/software/user/restricted/gurobi/11.0.2/licenses/gurobi.lic⁠
export HDF5_USE_FILE_LOCKING=FALSE
snakemake --cluster "sbatch -A {cluster.account} --mail-type ALL --mail-user {cluster.email}  -p {cluster.partition} -o {cluster.output} -e {cluster.error} -c {threads} --mem {resources.mem_mb} --time {cluster.walltime} --export=ALL,GRB_THREADS=8" --cluster-config config/config.cluster.yaml --jobs 20 --latency-wait 60 --configfile config/config.mi_rps.yaml --rerun-incomplete
