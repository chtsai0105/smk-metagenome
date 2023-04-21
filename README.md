[![Contributors][contributors-shield]][contributors-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]

[contributors-shield]: https://img.shields.io/github/contributors/chtsai0105/smk-metagenome
[contributors-url]: https://github.com/chtsai0105/smk-metagenome/graphs/contributors
[issues-shield]: https://img.shields.io/github/issues/chtsai0105/smk-metagenome
[issues-url]: https://github.com/chtsai0105/smk-metagenome/issues
[license-shield]: https://img.shields.io/github/license/chtsai0105/smk-metagenome?label=license
[license-url]: https://github.com/chtsai0105/smk-metagenome/blob/master/LICENSE.md

# Snakemake workflow for comparitive genomics
<img align="right" width="120" height="120" src="https://avatars.githubusercontent.com/u/33450111?s=200&v=4">
The main goal of this workflow is to perform taxonomic profiling and functional annotation for shotgun metagenome data.

The pair-ended shotgun sequencing data are first pre-processed by
[fastqc](https://github.com/s-andrews/FastQC)/[fastp](https://github.com/OpenGene/fastp)/[bfc](https://github.com/lh3/bfc) for
QC, trimming/adaptor removal and read error correction.
The processed reads are assembled by [metaSPAdes](https://github.com/ablab/spades) or
[megahit](https://github.com/voutcn/megahit) to obtain contigs.
The assembled contigs are further processed by [Autometa](https://github.com/KwanLab/Autometa) or
[metabat](https://bitbucket.org/berkeleylab/metabat)/[gtdbtk](https://github.com/Ecogenomics/GTDBTk) for taxonomic profiling.
The metabat outputs can also being sent to dbcan/kofamscan for functional annotation.

Currently the workflow is still under development. The results are not guaranteed correct.
Please verify all the results if you intend to publish them.

<br>

## Install the workflow management system [**snakemake**](https://snakemake.readthedocs.io/en/stable/index.html)
Snakemake is a Python based language and execution environment for GNU Make-like workflows.
The workflow is defined by specifying *rules* in a `snakefile`.
Rules further specify how to create sets of output files from sets of input files as well as the parameters and the command.

Snakemake automatically determines the dependencies between the rules by matching file names.
Please refer to [snakemake tutorial page](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) to see how to define the rules.

We have to make sure whether snakemake is installed in your system before we start.

### UCR hpcc
UCR hpcc have snakemake installed in environment module. If you're on UCR hpcc you can simply load snakemake module by the following command.
```
module load snakemake
```

### If you're not on UCR hpcc
If you're NOT on UCR hpcc and you don't have snakemake in the enviornment module, please follow the below steps to create a snakemake conda environment.
1. First, create an environment named **snakemake**.

    ```
    conda create -n snakemake
    ```

    After create the environment, activate it by:
    
    ```
    conda activate snakemake
    ```

2. Install the package **mamba**, which is a faster version of **conda**. 

    ```
    conda install -c conda-forge mamba
    ```
    
    After **mamba** being installed, you can later switch from `conda install [package]` to `mamba install [package]` to speed up the package installation.

3. Next, install the package **snakemake** through **mamba**.
    
    ```
    mamba install snakemake
    ```
    
4. Then you can execute the `snakemake --help` to show the snakemake helping page.

<br>

## Clone the workflow

Clone the repo to your computer.

Clone by the following command if you're using public key for github connection.

```
git clone --recurse-submodules git@github.com:chtsai0105/smk-metagenome.git
```

Or clone by https link.

```
git clone --recurse-submodules https://github.com/chtsai0105/smk-metagenome.git
```

Otherwise, clone the submodules as a second step by:
```
cd smk-metagenome
git submodule update --init
```

Next, go to the directory by `cd smk-metagenome`. It should contains the following files:

File    |Description
-|-
[snakefile](snakefile)  |The workflow entry. Define the targets for the workflow.
[config.yaml](config.yaml)  |Define the path for data and metadata.
[sample.csv](sample.csv)    |The metadata for samples. Define the names of the samples and the fastq files.
[run_snakemake.bash](run_snakemake.bash)    |The bash script for running the workflow.
[data/](data)   |The folder for the data and the workflow outputs.
[envs/](envs)   |The folder that contains the yaml config for conda environments.
[rules/](rules) |The folder that contains the rules/submodules of the workflow.
[slurm/](https://github.com/chtsai0105/snakemake_profile-slurm/tree/master) |The folder that contains the slurm profile for stajichlab partition@UCR hpcc.

<br>

## Configure the workflow

You can edit the `config.yaml` to setup the behavior of the workflow.

The key **metadata** refers to the `sample.csv`, which have all the details of the sample.
The other keys represent major steps in the workflow. You can switch on/off a particular step by changing the subkey `run` to True/False. You can also change the 
output path in some of the steps.

<br>

## Define the samples

You should properly defined your metadata, which is recorded in the `sample.csv`, before running the workflow.

There are 4 columns in this csv table - **sample**, **R1**, **R2** and **interleaved**.

The column **sample** defines the sample name. You can change it to names which are more distinguishable instead of accession numbers.
These names will also being used as the wildcard to define the workflow outputs.

The column **R1** and **R2** defines the fastq file names you put in the folder `data/fastq`.
If the sequencing data is in interleaved format. Please leave the **R1** and **R2** blanked and fill the interleaved file names to the **interleaved** column.

Please make sure they are identical to the fastq files you have otherwise the workflow may have trouble to input the files.
Please also confirm that the names in each column are unique.

<br>

## Run the workflow

After compiling the template and setup the paramters, the next step is to run the workflow.

Snakemake provides a dry-run feature to examine the workflow before truly running it. You should always test the workflow beforehand to make sure it execute as 
expected by the following command:

```
snakemake -np
```

After confirming all the steps. You can run the workflow by executing the script `run_snakemake.bash` or the following command:

```
snakemake -p --profile slurm --use-envmodules --use-conda --jobs 8 --max-threads 20
```

Currently I constrain the cpu usage of all parallel jobs to 20, you may change it as you like.
