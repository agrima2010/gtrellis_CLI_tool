# gtrellis_CLI_tool
Using gtrellis tool to create genome level trellis graph for visualization from the QC coverage output of epi2me-labs/wf-somatic-variation pipeline. This repository describes how to use an apptainer image to reproduce the analysis and use it as a CLI to generate QC plots in a closed HPC environment

# Reproducible Environment
This project was developed using an Apptainer (formerly Singularity) container built from a Conda environment. The container ensures reproducibility across HPC systems such as UPPMAX and other Linux clusters.

# Clone repository
```bash
git clone <repo>
cd gtrellis_CLI_tool
```
# Build the container
```bash
apptainer build r-gtrellis.sif Singularity.def
```
This creates an Apptainer image containing:

- R
- Bioconductor GTrellis
- GenomicRanges
- all other dependencies listed in environment.yml

# Usage
```bash
apptainer exec /path/to/file/r-gtrellis.sif Rscript plot_bed_tracks.R X_normal.regions.bed X_tumor.regions.bed
```
# Note
This project uses Apptainer to package the complete software stack, ensuring identical versions of R and Bioconductor packages across systems. This makes the workflow reproducible on HPC clusters where users may not have administrative privileges.


