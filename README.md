# gtrellis_CLI_tool
Using gtrellis tool to create genome level trellis graph for visualization from the QC coverage output of epi2me-labs/wf-somatic-variation pipeline. This repository describes how to use an apptainer image to reproduce the analysis and use it as a CLI to generate QC plots in a closed HPC environment

# Usage
apptainer exec /path/to/file/r-gtrellis.sif Rscript plot_bed_tracks.R X_normal.regions.bed X_tumor.regions.bed


