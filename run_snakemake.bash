#!/bin/bash

snakemake -p --profile slurm --use-envmodules --use-conda --jobs 8 --max-threads 20
