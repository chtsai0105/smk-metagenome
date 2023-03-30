import os
import pandas as pd
import re


class df_prebuild(object):
    def __init__(self, df):
        df = df.copy()
        self.__check(df)
        df = self.__transform(df)
        df = self.__add_content(df)
        
        self.df = df
    
    def __check(self, df):
        for idx, v in df.iterrows():
            if bool(v['R1']) ^ bool(v['R2']):
                raise ValueError("One of the paired-end reads is missing (R1 or R2).")
            if not bool(v['R1']) ^ bool(v['interleaved']):
                raise ValueError("R1/R2 and interleaved should be mutually exclusive.")

    def __transform(self, df):
        return df.melt(id_vars=['sample', 'interleaved'], var_name='read', value_vars=['R1', 'R2'], value_name='origin')
    
    def __add_content(self, df):
        '''
        Match extensions for .fastq, .fq, .fastq.gz, .fq.gz
        '''
        temp = df.query('interleaved != ""')
        df.loc[df['interleaved'] != "", 'origin'] = temp['sample'] + '_deinterleaved_' + temp['read'] + '.fastq.gz'
        df['ext'] = df['origin'].apply(lambda x: re.search("\.f(ast)?q($|\.gz$)", x).group())
        df['base'] = df['sample'] + '_' + df['read']
        df['renamed_fastq'] = df['base'] + df['ext']
        return df
    
    def samples(self):
        return self.df['sample'].drop_duplicates().tolist()
    
    def deinterleaved_samples(self):
        df = self.df.query('interleaved != ""')
        return df['origin'].tolist()
    
    def renamed_samples(self):
        df = self.df
        return df.query('origin != renamed_fastq')['renamed_fastq'].tolist()
        
    def unified_samples(self, sample, path, read=None, column="base", ext=""):
        '''
        Return unified sample names
        Parameter:
            sample: sample name. Take wildcards.sample in snakemake workflow.
            path: the path of the parent folder of the files.
            read: direction of the read. Should be R1 or R2. Not specifying will return a list containing both R1 and R2.
            column: which column in the df.df you want to return. By default it returns the basename column.
            ext: a custom extension you want to added after the file. Should be used with column="base".
        '''
        if ext and column != "base":
            raise ValueError("When specifying the variable ext, the variable column should be \"base\"")
        if read:
            return [os.path.join(path, x + '{}'.format(ext)) for x in self.df.query('sample == @sample').query('read == @read')[column].tolist()]
        else:
            return [os.path.join(path, x + '{}'.format(ext)) for x in self.df.query('sample == @sample')[column].drop_duplicates().tolist()]

    
############### Configuration ##############
configfile: "config.yaml"
data = pd.read_csv(config['Metadata'], keep_default_na=False, na_values=['_'], comment="#")
data = df_prebuild(data)

FASTQ = config['raw_data']['fastq']
FASTQC_OUTPUT = config['fastqc']['output']
FASTQ_TRIMMED = config['trimming']['output']
KRAKEN2_OUTPUT = config['kraken2']['output']
ASSEMBLY_OUTPUT = config['assembly']['output']
FILTERED_CONTIGS = config['assembly']['filtered_contigs']
MAPPING_OUTPUT = config['align_to_assembly']['output']
AUTOMETA_OUTPUT = config['autometa']['output']
METABAT_OUTPUT = config['metabat']['output']
GTDBTK_OUTPUT = config['gtdbtk']['output']
FUNC_ANNO_OUTPUT = config['functional_annotation']['output']
############### Input settings #############
input_list = list()

### Deinterleave and rename fastq
input_list.extend(["{dir}/{sample}".format(dir=FASTQ, sample=sample) for sample in data.deinterleaved_samples()])
input_list.extend(["{dir}/{sample}".format(dir=FASTQ, sample=sample) for sample in data.renamed_samples()])

### FastQC
if config['fastqc']['run']:
    input_list.extend(["{dir}/pre_trim/{sample}_R1_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/pre_trim/{sample}_R1_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/pre_trim/{sample}_R2_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/pre_trim/{sample}_R2_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])

### Trimmomatic and post-trim fastqc
if config['trimming']['run']:
    input_list.extend(["{dir}/{sample}_R1.fastq.gz".format(dir=FASTQ_TRIMMED, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/{sample}_R2.fastq.gz".format(dir=FASTQ_TRIMMED, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/post_trim/{sample}_R1_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/post_trim/{sample}_R1_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/post_trim/{sample}_R2_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/post_trim/{sample}_R2_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in data.samples()])

### Kraken2
if config['kraken2']['run']:
    include: "rules/kraken2.smk"
    input_list.extend(["{dir}/{sample}_taxon.csv".format(dir=KRAKEN2_OUTPUT, sample=sample) for sample in data.samples()])

### Assembly

if config['assembly']['run']:
    include: "rules/assembly.smk"
    input_list.extend(["{dir}/{sample}_contigs.fasta".format(dir=FILTERED_CONTIGS, sample=sample) for sample in data.samples()])
    input_list.extend(["{dir}/{sample}_filtered.fasta".format(dir=FILTERED_CONTIGS, sample=sample) for sample in data.samples()])  # filtered_fasta

### Post-assembly alignment
if config['align_to_assembly']['run']:
    include: "rules/post_checkup.smk"
    if config['align_to_assembly']['tools'] == "bbmap":
        input_list.extend(["{dir}/{sample}_refstats.tsv".format(dir=MAPPING_OUTPUT, sample=sample) for sample in data.samples()])
    else:
        input_list.extend(["{dir}/{sample}.stats".format(dir=MAPPING_OUTPUT, sample=sample) for sample in data.samples()])

### euk_detection.smk
# if config['run_euk_detection']:
#     include: "rules/euk_detection.smk"
#     input_list.extend(["{dir}/{sample}/bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in data.samples()])
#     input_list.extend(["{dir}/{sample}/euk_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in data.samples()])
#     input_list.extend(["{dir}/{sample}/prok_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in data.samples()])

### autometa.smk
if config['autometa']['run']:
    include: "rules/autometa.smk"

    input_list.extend(["{dir}/{sample}/intermediates/coverage.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in data.samples()])  # cov_tab
    input_list.extend(["{dir}/{sample}/intermediates/blastp.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in data.samples()])    # blastp
    input_list.extend(["{dir}/{sample}/intermediates/taxonomy/taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in data.samples()])   # taxonomy
    
    for kingdom in config['autometa']['binning_target']:
        # input_list.extend(["{dir}/{sample}/intermediates/{kingdom}.markers.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()])  # autometa_markers
        # input_list.extend(["{dir}/{sample}/{kingdom}_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()])  # binning_output
        # input_list.extend(["{dir}/{sample}/{kingdom}_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()]) # main_output
        if config['autometa']['unclustered_recruitment']:
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()])    # metabin_stats
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_features.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()]) # metabin_taxonomy
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()]) # metabin
        input_list.extend(["{dir}/{sample}/{kingdom}_metabin_stats.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()])    # metabin_stats
        input_list.extend(["{dir}/{sample}/{kingdom}_metabin_taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()]) # metabin_taxonomy
        input_list.extend(["{dir}/{sample}/{kingdom}_metabins".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in data.samples()]) # metabin

if config['metabat']['run']:
    include: "rules/metabat.smk"
    input_list.extend(["{dir}/{sample}/bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in data.samples()])

if config['gtdbtk']['run']:
    include: "rules/gtdbtk.smk"
    input_list.extend(["{dir}/{sample}".format(dir=GTDBTK_OUTPUT, sample=sample) for sample in data.samples()])

if config['functional_annotation']['run']:
    include: "rules/functional_annotation.smk"
    if config['functional_annotation']['clustering']:
        input_list.extend(["{dir}/intermediate/{sample}_alignment_readcounts.csv".format(dir=FUNC_ANNO_OUTPUT, sample=sample) for sample in data.samples()])
    if "dbcan" in config['functional_annotation']['tools']:
        input_list.extend(["{dir}/dbcan/{sample}/overview.txt".format(dir=FUNC_ANNO_OUTPUT, sample=sample) for sample in data.samples()])
    if "kofamscan" in config['functional_annotation']['tools']:
        input_list.extend(["{dir}/kofamscan/{sample}.txt".format(dir=FUNC_ANNO_OUTPUT, sample=sample) for sample in data.samples()])
    


############### Rules ######################
localrules: rename_input

ruleorder: deinterleave > rename_input

wildcard_constraints:
        sample = "[^/]+",                # Regex for all characters except /
        R = "R(1|2)",
        ext = "f(ast)?q($|\.gz$)"


rule all:
    input:
        input_list

rule deinterleave:
    input:
        lambda w: data.unified_samples(w.sample, FASTQ, column='interleaved')
    output:
        R1 = "{dir}/{{sample}}_deinterleaved_R1.fastq.gz".format(dir=FASTQ),
        R2 = "{dir}/{{sample}}_deinterleaved_R2.fastq.gz".format(dir=FASTQ)
    wildcard_constraints:
        sample = "[^/]+[^(_deinterleaved)]"
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in={input} out1={output.R1} out2={output.R2}
        """

rule rename_input:
    input:
        lambda w: data.unified_samples(w.sample, FASTQ, read=w.R, column='origin')
    output:
        "{dir}/{{sample}}_{{R}}.{{ext}}".format(dir=FASTQ)
    shell:
        """
        ln -sr {input} {output}
        """

rule fastqc_pre:
    input:
        lambda w: data.unified_samples(w.sample, FASTQ, read=w.R, column='renamed_fastq')
    output:
        expand("{dir}/pre_trim/{{sample}}_{{R}}_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"])
    params:
        dirname = "{dir}/pre_trim".format(dir=FASTQC_OUTPUT)
    threads: 4
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.dirname} {input}
        """

rule trimmomatic:
    input:
        R1 = lambda w: data.unified_samples(w.sample, FASTQ, read='R1', column='renamed_fastq'),
        R2 = lambda w: data.unified_samples(w.sample, FASTQ, read='R2', column='renamed_fastq')
    output:
        R1_paired = "{dir}/{{sample}}_R1.fastq.gz".format(dir=FASTQ_TRIMMED),
        R1_unpaired = "{dir}/{{sample}}_R1_unpaired.fastq.gz".format(dir=FASTQ_TRIMMED),
        R2_paired = "{dir}/{{sample}}_R2.fastq.gz".format(dir=FASTQ_TRIMMED),
        R2_unpaired = "{dir}/{{sample}}_R2_unpaired.fastq.gz".format(dir=FASTQ_TRIMMED)
    threads: 4
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        trimmomatic PE -threads {threads} {input.R1} {input.R2} \
        {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
        """

use rule fastqc_pre as fastqc_post with:
    input:
        lambda w: rules.trimmomatic.output.R1_paired if w.R == 'R1'
            else rules.trimmomatic.output.R2_paired
    output:
        expand("{dir}/post_trim/{{sample}}_{{R}}_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"])
    params:
        dirname = "{dir}/post_trim".format(dir=FASTQC_OUTPUT)
