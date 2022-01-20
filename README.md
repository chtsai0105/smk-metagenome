## Install the workflow manage system [**snakemake**](https://snakemake.readthedocs.io/en/stable/index.html)
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
    
4. Then you can execute the `snakemake --help` to show the snakemake helping page. Snakemake is a Python based language and execution environment for GNU Make-like workflows. The workflow is defined by specifying rules in a `snakefile`. Rules further specify how to create sets of output files from sets of input files as well as the parameters and the command. Snakemake automatically determines the dependencies between the rules by matching file names. Please refer to [snakemake tutorial page](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html) to see how to define the rules.

## Clone the metagenome workflow

Clone the repo to your computer.

Clone by the following command if you're using public key for github connection.

```
git clone git@github.com:chtsai0105/metagenome-snakemake.git
```

Otherwise, clone by https link.

```
git clone https://github.com/chtsai0105/metagenome-snakemake.git
```

Next, go to the directory by `cd metagenome-snakemake`. It should contains the following files:

File                    |Description
------------------------|---------------------------------
`snakefile`             |Define the rules for the workflow.
`config.yaml`           |Define the path for data and metadata.
`sample.csv`            |The metadata for samples. Define the names of the samples and the fastq files.
`data`                  |The folder for the data and the workflow outputs.
`run_snakemake.bash`    |Execute the workflow

## Define the path

You can edit the `config.yaml` to setup the path for data and metadata of the samples. By default, the metadata is defined in the file `sample.csv` and all the data should be stored in the folder `data`.

## Define the samples

You should properly defined your metadata, which is record in the `sample.csv`, before running the workflow. There are 2 columns in this csv table - **sample** and **fastq**. The column **sample** defined the sample name. You can change it to names which are more distinguishable instead of a accession number. The column **fastq** defined the fastq file names. Please make sure they are identical to the fastq files you have otherwise the workflow may have trouble to input the files. Please also confirm that the names in each column are unique.

## Run the workflow

After compiling the template and setup the paramters, the next step is to run the workflow.

Snakemake provide a dry-run feature to examine the workflow before truly running it. You should always use the following command to test the workflow beforehand:

```
snakemake -np
```

After confirming all the steps. You can run the workflow by executing the script `run_snakemake.bash` or the following command:

```
snakemake -p --profile slurm --use-conda
```
