import os
import pandas as pd
import re

def extract_ext(path_str):
    '''
    Match extensions for .fastq, .fq, .fastq.gz, .fq.gz
    '''
    m = re.search("\.f(ast)?q($|\.gz$)", path_str)
    ext = m.group()
    return ext

def fastq_input(wildcards, p, r):
    if r not in ('R1', 'R2', 'interleaved'):
        print("Error, only accept R1, R2 or interleaved")
    return os.path.join(p, sample_df.loc[sample_df['sample'] == wildcards.sample, r].item())

def fastqc_input(wildcards, p):
    r = re.compile("{}.*".format(wildcards.sample))
    f = list(filter(r.match, sample_df[['R1', 'R2']].melt()['value']))[0]
    return os.path.join(p, f)

def trimmed_fastq_input(wildcards, p, r):
    return os.path.join(p, sample_df.loc[sample_df['sample'] == wildcards.sample, 'sample'].item()) + '_{}.fastq.gz'.format(r)

    
############### Configuration ##############
configfile: "config.yaml"
sample_df = pd.read_csv(config['Metadata'], keep_default_na=False, na_values=['_'], comment="#")
for idx, v in sample_df.iterrows():
    if v['R1'] != "" and v['R2'] != "" and v['interleaved'] == "":
        sample_df.loc[idx, 'ext'] = extract_ext(v['R1'])
    elif v['R1'] == "" and v['R2'] == "" and v['interleaved'] != "":
        sample_df.loc[idx, 'ext'] = '.fastq.gz'
        sample_df.loc[idx, 'R1'] = sample_df.loc[idx, 'sample'] + '_deinterleaved_R1' + sample_df.loc[idx, 'ext']
        sample_df.loc[idx, 'R2'] = sample_df.loc[idx, 'sample'] + '_deinterleaved_R2' + sample_df.loc[idx, 'ext']


############### Input settings #############
input_list = list()

FASTQ = config['raw_data']['fastq']
input_list.extend(["{dir}/{fastq}".format(dir=FASTQ, fastq=fastq) for fastq in sample_df['R1']])
input_list.extend(["{dir}/{fastq}".format(dir=FASTQ, fastq=fastq) for fastq in sample_df['R2']])

FASTQC_OUTPUT = config['fastqc']['output']
if config['fastqc']['run_pretrim_qc']:
    R1_base = sample_df.apply(lambda row : row['R1'].replace(row['ext'], ''), axis=1).tolist()
    R2_base = sample_df.apply(lambda row : row['R2'].replace(row['ext'], ''), axis=1).tolist()
    input_list.extend(["{dir}/pre_trim/{fastq}_fastqc.html".format(dir=FASTQC_OUTPUT, fastq=fastq) for fastq in R1_base + R2_base])
    input_list.extend(["{dir}/pre_trim/{fastq}_fastqc.zip".format(dir=FASTQC_OUTPUT, fastq=fastq) for fastq in R1_base + R2_base])

### Trimmomatic and post-trim fastqc
FASTQ_TRIMMED = config['trimming']['output']
FASTQC_OUTPUT = config['fastqc']['output']
if config['trimming']['run_trimmomatic']:
    input_list.extend(["{dir}/{fastq}_R1.fastq.gz".format(dir=FASTQ_TRIMMED, fastq=fastq) for fastq in sample_df['sample']])
    input_list.extend(["{dir}/{fastq}_R2.fastq.gz".format(dir=FASTQ_TRIMMED, fastq=fastq) for fastq in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R1_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R1_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R2_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R2_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])

KRAKEN2_OUTPUT = config['kraken2']['output']
if config['kraken2']['run_kraken2']:
    include: "rules/kraken2.smk"
    input_list.extend(["{dir}/{sample}_taxon.csv".format(dir=KRAKEN2_OUTPUT, sample=sample) for sample in sample_df['sample']])

### Assembly
ASSEMBLY_OUTPUT = config['assembly']['output']
FILTERED_CONTIGS = config['assembly']['filtered_contigs']
if config['assembly']['run_assembly']:
    include: "rules/assembly.smk"
    input_list.extend(["{dir}/{sample}_contigs.fasta".format(dir=FILTERED_CONTIGS, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/{sample}_filtered.fasta".format(dir=FILTERED_CONTIGS, sample=sample) for sample in sample_df['sample']])  # filtered_fasta

# if config['postassembled_alignment']['run_alignment']:
#     input_list.extend(["{dir}/{sample}_covstats.tsv".format(dir=MAPPING_OUTPUT, sample=sample) for sample in sample_df['sample']])

### euk_detection.smk
# if config['run_euk_detection']:
#     include: "rules/euk_detection.smk"
#     input_list.extend(["{dir}/{sample}/bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])
#     input_list.extend(["{dir}/{sample}/euk_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])
#     input_list.extend(["{dir}/{sample}/prok_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])

### autometa.smk
AUTOMETA_OUTPUT = config['autometa']['output']
if config['autometa']['run_autometa']:
    include: "rules/autometa.smk"

    input_list.extend(["{dir}/{sample}/intermediates/coverage.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])  # cov_tab
    input_list.extend(["{dir}/{sample}/intermediates/blastp.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])    # blastp
    input_list.extend(["{dir}/{sample}/intermediates/taxonomy/taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])   # taxonomy
    
    for kingdom in config['autometa']['binning_target']:
        # input_list.extend(["{dir}/{sample}/intermediates/{kingdom}.markers.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])  # autometa_markers
        # input_list.extend(["{dir}/{sample}/{kingdom}_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])  # binning_output
        # input_list.extend(["{dir}/{sample}/{kingdom}_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # main_output
        if config['autometa']['unclustered_recruitment']:
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])    # metabin_stats
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_features.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin_taxonomy
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin
        input_list.extend(["{dir}/{sample}/{kingdom}_metabin_stats.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])    # metabin_stats
        input_list.extend(["{dir}/{sample}/{kingdom}_metabin_taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin_taxonomy
        input_list.extend(["{dir}/{sample}/{kingdom}_metabins".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin

MAPPING_OUTPUT = config['post_checkup']['output']
if config['post_checkup']['align_MAGs']:
    include: "rules/post_checkup.smk"
    if config['post_checkup']['tools'] == 'bbmap':
        input_list.extend(["{dir}/{sample}_refstats.tsv".format(dir=MAPPING_OUTPUT, sample=sample) for sample in sample_df['sample']])
    else:
        input_list.extend(["{dir}/{sample}.stats".format(dir=MAPPING_OUTPUT, sample=sample) for sample in sample_df['sample']])


print(input_list)
print(sample_df)
############### Rules ######################
wildcard_constraints:
        sample = "[^/]+"                # Regex for all characters except /

rule all:
    input:
        input_list

rule deinterleave:
    input:
        lambda wildcards: fastq_input(wildcards, FASTQ, 'interleaved')
    output:
        R1 = "{dir}/{{sample}}_deinterleaved_R1.fastq.gz".format(dir=FASTQ),
        R2 = "{dir}/{{sample}}_deinterleaved_R2.fastq.gz".format(dir=FASTQ)
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in={input} out1={output.R1} out2={output.R2}
        """

rule fastqc_pre:
    input:
        lambda wildcards: fastqc_input(wildcards, FASTQ)
    output:
        expand("{dir}/pre_trim/{{sample}}_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"])
    params:
        dirname = "{dir}/pre_trim".format(dir=FASTQC_OUTPUT)
    threads: 4
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.dirname} {input}
        """

##### Trimmomatic and post-trim FastQC #####
rule trimmomatic:
    input:
        R1 = lambda wildcards: fastq_input(wildcards, FASTQ, 'R1'),
        R2 = lambda wildcards: fastq_input(wildcards, FASTQ, 'R2')
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
        trimmomatic PE -threads {threads} {input.R1} {input.R2} {output.R1_paired} {output.R1_unpaired} {output.R2_paired} {output.R2_unpaired} \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
        """

use rule fastqc_pre as fastqc_post with:
    input:
        lambda wildcards: rules.trimmomatic.output.R1_paired if wildcards.read == 'R1' else rules.trimmomatic.output.R2_paired
    output:
        expand("{dir}/post_trim/{{sample}}_{{read}}_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"])
    params:
        dirname = "{dir}/post_trim".format(dir=FASTQC_OUTPUT)
############################################
