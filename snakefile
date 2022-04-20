import os
import pandas as pd


def extract_ext(path_str):
    string, ext = os.path.splitext(path_str)
    if ext == ".gz":
        ext = os.path.splitext(string)[-1] + ext
    return ext  # Avoid the "." at the beginning of the ext string
    

configfile: "config.yaml"
sample_df = pd.read_csv(config['Metadata'])
sample_df['ext'] = sample_df['fastq'].apply(extract_ext)
sample_df['fastq_renamed'] = sample_df['sample'] + sample_df['ext']

############# Program settings #############

run_trimmomatic = False
include: "includes/euk_detection.smk"
include: "includes/autometa.smk"

############################################

############### Path settings ##############
DATA = config['Data']
FASTQ = os.path.join(DATA, 'fastq')
FASTQ_RENAMED = os.path.join(DATA, 'fastq_renamed')
if run_trimmomatic:
    FASTQ_TRIMMED = os.path.join(DATA, 'fastq_trimmed')
else:
    FASTQ_TRIMMED = FASTQ_RENAMED
FASTQC_OUTPUT = os.path.join(DATA, 'fastqc')

ASSEMBLY_OUTPUT = os.path.join(DATA, 'spades')
MAPPING_OUTPUT = os.path.join(DATA, 'bowtie2')
BINNING_OUTPUT = os.path.join(DATA, 'autometa')
EUK_BINNING_OUTPUT = os.path.join(DATA, 'metabat_euk')

############################################

localrules: rename_input
wildcard_constraints:
        ext = "f(ast)?q($|\.gz$)",      # Regex for fastq, fq, fastq.gz and fq.gz as extension
        sample = "[^/]+"                # Regex for all characters except /

rule all:
    input:
        fastq_renamed =     expand("{dir}/{fastq}", dir=FASTQ_RENAMED, fastq=sample_df['fastq_renamed']),
        fastqc_pre_trim =   expand("{dir}/pre_trim/{sample}_fastqc.{ext}", dir=FASTQC_OUTPUT, sample=sample_df['sample'], ext=["html", "zip"]),
        # fastq_trimmed =     expand("{dir}/{fastq}", dir=FASTQ_TRIMMED, fastq=sample_df['fastq_renamed']),
        # fastqc_post_trim =  expand("{dir}/post_trim/{sample}_fastqc.{ext}", dir=FASTQC_OUTPUT, sample=sample_df['sample'], ext=["html", "zip"]),
        assembly =          expand("{dir}/{sample}/scaffolds.fasta", dir=ASSEMBLY_OUTPUT, sample=sample_df['sample']),
        mapping_output =    expand("{dir}/{sample}.{ext}", dir=MAPPING_OUTPUT, sample=sample_df['sample'], ext=["bam", "bam.bai"]),
        ### euk_detection.smk 
        # euk_bin =           expand("{dir}/{sample}/euk_bin", dir=EUK_BINNING_OUTPUT, sample=sample_df['sample']),
        #####
        ### autometa.smk
        filtered_fasta =    expand("{dir}/{sample}/intermediates/filtered.fasta", dir=BINNING_OUTPUT, sample=sample_df['sample']),
        # cov_tab =           expand("{dir}/{sample}/coverage.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample']),
        # autometa_markers =  expand("{dir}/{sample}/{kingdom}.markers.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=["bacteria", "archaea"]),
        # blastp =            expand("{dir}/{sample}/blastp.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample']),
        # kingdom_fasta =     expand("{dir}/{sample}/{kingdom}.fasta", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=["bacteria", "archaea"]),
        taxonomy =          expand("{dir}/{sample}/intermediates/taxonomy/taxonomy.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample']),
        # binning_output =    expand("{dir}/{sample}/{kingdom}_binning.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea']),
        # main_output =       expand("{dir}/{sample}/{kingdom}_main.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea']),
        # recruit_binning =   expand("{dir}/{sample}/{kingdom}_recruitment_binning.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea']),
        # recruit_features =  expand("{dir}/{sample}/{kingdom}_recruitment_features.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea']),
        # recruit_main =      expand("{dir}/{sample}/{kingdom}_recruitment_main.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea']),
        metabin_stats = expand("{dir}/{sample}/{kingdom}_metabin_stats.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea']),
        metabin_taxonomy = expand("{dir}/{sample}/{kingdom}_metabin_taxonomy.tsv", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea']),
        metabins = expand("{dir}/{sample}/{kingdom}_metabins", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['bacteria', 'archaea'])
        ###

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
if run_trimmomatic:
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
        lambda wildcards: os.path.join(FASTQ_TRIMMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
    output:
        assembly = "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT)
    params:
        dirname = directory("{dir}/{{sample}}".format(dir=ASSEMBLY_OUTPUT))
    threads: 12
    resources:
        time=config['time']['LV2'],
        mem_mb=config['mem']['LV3']
    conda:
        "envs/assembler.yaml"
    shell:
        """
        spades.py --meta -o {params.dirname} --12 {input} -t {threads}
        """


rule bowtie2_mapping:
    input:
        assembly = rules.metaspades.output.assembly,
        fastq = lambda wildcards: os.path.join(FASTQ_TRIMMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
    output:
        idx = temp(expand("{dir}/{{sample}}.{ext}.bt2", dir=MAPPING_OUTPUT, ext=["1", "2", "3", "4", "rev.1", "rev.2"])),
        bam = "{dir}/{{sample}}.bam".format(dir=MAPPING_OUTPUT),
        bai = "{dir}/{{sample}}.bam.bai".format(dir=MAPPING_OUTPUT)
    params:
        idx = "{dir}/{{sample}}".format(dir=MAPPING_OUTPUT)
    threads: 4
    resources:
        time=config['time']['LV1'],
        mem_mb=config['mem']['LV2']
    conda:
        "envs/autometa.yaml"
    shell:
        """
        bowtie2-build {input.assembly} {params.idx}
        bowtie2 -p {threads} -x {params.idx} --interleaved {input.fastq} | samtools view -@ {threads} -Sbhu - | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam} {output.bai}
        """
