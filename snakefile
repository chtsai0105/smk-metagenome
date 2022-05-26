import os
import pandas as pd


def extract_ext(path_str):
    ext = Path(path_str).suffixes
    return "".join(ext)
    

configfile: "config.yaml"
sample_df = pd.read_csv(config['Metadata'])
sample_df['ext'] = sample_df['fastq'].apply(extract_ext)
sample_df['fastq_renamed'] = sample_df['sample'] + sample_df['ext']

############# Program settings #############


include: "rules/euk_detection.smk"
include: "rules/autometa.smk"

############################################

############### Path settings ##############
FASTQ = config['Path']['fastq']
FASTQ_RENAMED = config['Path']['fastq_renamed']

if config['run_trimmomatic']:
    FASTQ_TRIMMED = config['Path']['fastq_trimmed']
else:
    FASTQ_TRIMMED = FASTQ_RENAMED
FASTQC_OUTPUT = config['Path']['fastqc_output']

ASSEMBLY_OUTPUT = config['Path']['assembly_output']
MAPPING_OUTPUT = config['Path']['mapping_output']
AUTOMETA_OUTPUT = config['Path']['autometa_output']
METABAT_OUTPUT = config['Path']['metabat_output']

############################################

############### Input settings #############

input_list = list()

input_list.extend(["{dir}/{fastq}".format(dir=FASTQ_RENAMED, fastq=fastq) for fastq in sample_df['fastq_renamed']])
input_list.extend(["{dir}/pre_trim/{sample}_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
input_list.extend(["{dir}/pre_trim/{sample}_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])

### Trimmomatic and post-trim fastqc
if config['run_trimmomatic']:
    input_list.extend(["{dir}/{fastq}".format(dir=FASTQ_TRIMMED, fastq=fastq) for fastq in sample_df['fastq_renamed']])
    input_list.extend(["{dir}/post_trim/{sample}_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])

### Assembly
input_list.extend(["{dir}/{sample}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT, sample=sample) for sample in sample_df['sample']])

### euk_detection.smk
if config['run_euk_detection']:
    input_list.extend(["{dir}/{sample}.bam".format(dir=MAPPING_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/{sample}.bam.bai".format(dir=MAPPING_OUTPUT, sample=sample) for sample in sample_df['sample']])
    # input_list.extend(["{dir}/{sample}/bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/{sample}/euk_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/{sample}/prok_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])

### autometa.smk
if config['autometa']['run_autometa']:
    input_list.extend(["{dir}/{sample}/intermediates/filtered.fasta".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])  # filtered_fasta
    input_list.extend(["{dir}/{sample}/intermediates/coverage.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])  # cov_tab
    input_list.extend(["{dir}/{sample}/intermediates/blastp.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])    # blastp
    input_list.extend(["{dir}/{sample}/intermediates/taxonomy/taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])   # taxonomy
    
    for kingdom in config['autometa']['binning_target']:
        input_list.extend(["{dir}/{sample}/intermediates/{kingdom}.markers.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])  # autometa_markers
        input_list.extend(["{dir}/{sample}/{kingdom}_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])  # binning_output
        input_list.extend(["{dir}/{sample}/{kingdom}_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # main_output
        if config['autometa']['unclustered_recruitment']:
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])    # metabin_stats
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_features.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin_taxonomy
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin
        input_list.extend(["{dir}/{sample}/{kingdom}_metabin_stats.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])    # metabin_stats
        input_list.extend(["{dir}/{sample}/{kingdom}_metabin_taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin_taxonomy
        input_list.extend(["{dir}/{sample}/{kingdom}_metabins".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin

############################################

localrules: rename_input
wildcard_constraints:
        ext = "f(ast)?q($|\.gz$)",      # Regex for fastq, fq, fastq.gz and fq.gz as extension
        sample = "[^/]+"                # Regex for all characters except /

rule all:
    input:
        input_list

rule rename_input:
    input:
        lambda wildcards: os.path.abspath(os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq'].item()))
    output:
        "{dir}/{{sample}}.{{ext}}".format(dir=FASTQ_RENAMED)
    shell:
        """
        ln -s {input} {output}
        """

rule fastqc_pre:
    input:
        lambda wildcards: os.path.join(FASTQ_RENAMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
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
        rules.rename_input.output
    output:
        "{dir}/{{sample}}.{{ext}}".format(dir=FASTQ_TRIMMED)
    threads: 4
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        trimmomatic SE -threads {threads} {input} {output} \
        ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15
        """

use rule fastqc_pre as fastqc_post with:
    input:
        lambda wildcards: os.path.join(FASTQ_TRIMMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
    output:
        expand("{dir}/post_trim/{{sample}}_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"])
    params:
        dirname = "{dir}/post_trim".format(dir=FASTQC_OUTPUT)

############################################

rule metaspades:
    input:
        lambda wildcards: os.path.join(FASTQ_TRIMMED if config['run_trimmomatic'] else FASTQ_RENAMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
    output:
        assembly = "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT)
    params:
        dirname = directory("{dir}/{{sample}}".format(dir=ASSEMBLY_OUTPUT))
    threads: 12
    resources:
        time="14-00:00:00",
        mem_mb=lambda wildcards, input, attempt: (input.size // 1000000) * attempt * 100
    conda:
        "envs/assembler.yaml"
    shell:
        """
        spades.py --meta -o {params.dirname} --12 {input} -t {threads}
        """

rule filter_contig_length:
    input:
        rules.metaspades.output.assembly
    output:
        "{dir}/{{sample}}_filtered.fasta".format(dir=MAPPING_OUTPUT)
    conda:
        "envs/assembler.yaml"
    shell:
        """
        reformat.sh in={input} out={output} minlength=3000
        """

rule bowtie2_mapping:
    input:
        assembly = rules.filter_contig_length.output,
        fastq = lambda wildcards: os.path.join(FASTQ_TRIMMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
    output:
        idx = temp(expand("{dir}/{{sample}}.{ext}.bt2", dir=MAPPING_OUTPUT, ext=["1", "2", "3", "4", "rev.1", "rev.2"])),
        bam = "{dir}/{{sample}}.bam".format(dir=MAPPING_OUTPUT),
        bai = "{dir}/{{sample}}.bam.bai".format(dir=MAPPING_OUTPUT)
    params:
        idx = "{dir}/{{sample}}".format(dir=MAPPING_OUTPUT)
    threads: 8
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: (input.size // 1000000) * attempt * 20
    conda:
        "envs/autometa.yaml"
    shell:
        """
        bowtie2-build {input.assembly} {params.idx}
        bowtie2 -p {threads} -x {params.idx} --interleaved {input.fastq} | samtools view -@ {threads} -Sbhu - | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam} {output.bai}
        """
