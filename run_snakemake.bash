#!/bin/bash

snakemake -p --profile slurm --use-conda --jobs 6 --max-threads 20 --conda-prefix $HOME/.conda/snakemake
