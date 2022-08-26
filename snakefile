import os
import pandas as pd
import re

def extract_ext(path_str):
    m = re.search("\.f(ast)?q($|\.gz$)", path_str)
    ext = m.group()
    return ext
    
############### Configuration ##############
configfile: "config.yaml"
sample_df = pd.read_csv(config['Metadata'], keep_default_na=False, na_values=['_'], comment="#")
for idx, v in sample_df.iterrows():
    if v['R1'] != "" and v['R2'] != "" and v['interleaved'] == "":
        sample_df.loc[idx, 'ext'] = extract_ext(v['R1'])
    elif v['R1'] == "" and v['R2'] == "" and v['interleaved'] != "":
        sample_df.loc[idx, 'ext'] = extract_ext(v['interleaved'])
        sample_df.loc[idx, 'R1'] = sample_df.loc[idx, 'sample'] + '_R1' + sample_df.loc[idx, 'ext']
        sample_df.loc[idx, 'R2'] = sample_df.loc[idx, 'sample'] + '_R2' + sample_df.loc[idx, 'ext']


############### Path settings ##############
FASTQ = config['Path']['fastq']

FASTQ_TRIMMED = config['Path']['fastq_trimmed']
FASTQC_OUTPUT = config['Path']['fastqc_output']

SPADES_OUTPUT = config['Path']['spades_output']
MEGAHIT_OUTPUT = config['Path']['megahit_output']
FILTERED_CONTIGS = config['Path']['filtered_contigs']
MAPPING_OUTPUT = config['Path']['mapping_output']
AUTOMETA_OUTPUT = config['Path']['autometa_output']
METABAT_OUTPUT = config['Path']['metabat_output']


############### Input settings #############
input_list = list()
input_list.extend(["{dir}/{fastq}".format(dir=FASTQ, fastq=fastq) for fastq in sample_df['R1']])
input_list.extend(["{dir}/{fastq}".format(dir=FASTQ, fastq=fastq) for fastq in sample_df['R2']])
# input_list.extend(["{dir}/{fastq}".format(dir=FASTQ_RENAMED, fastq=fastq) for fastq in sample_df['fastq_renamed']])
input_list.extend(["{dir}/pre_trim/{sample}_R1_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
input_list.extend(["{dir}/pre_trim/{sample}_R1_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
input_list.extend(["{dir}/pre_trim/{sample}_R2_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
input_list.extend(["{dir}/pre_trim/{sample}_R2_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])

### Trimmomatic and post-trim fastqc
if config['run_trimmomatic']:
    # input_list.extend(["{dir}/{fastq}_R1.fastq.gz".format(dir=FASTQ_TRIMMED, fastq=fastq) for fastq in sample_df['sample']])
    # input_list.extend(["{dir}/{fastq}_R2.fastq.gz".format(dir=FASTQ_TRIMMED, fastq=fastq) for fastq in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R1_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R1_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R2_fastqc.html".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/post_trim/{sample}_R2_fastqc.zip".format(dir=FASTQC_OUTPUT, sample=sample) for sample in sample_df['sample']])

### Assembly
input_list.extend(["{dir}/{sample}_contigs.fasta".format(dir=FILTERED_CONTIGS, sample=sample) for sample in sample_df['sample']])
input_list.extend(["{dir}/{sample}_filtered.fasta".format(dir=FILTERED_CONTIGS, sample=sample) for sample in sample_df['sample']])  # filtered_fasta

# if config['align_against_scaffold']:
#     input_list.extend(["{dir}/{sample}_covstats.tsv".format(dir=MAPPING_OUTPUT, sample=sample) for sample in sample_df['sample']])

### euk_detection.smk
if config['run_euk_detection']:
    include: "rules/euk_detection.smk"
    # input_list.extend(["{dir}/{sample}/bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/{sample}/euk_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])
    input_list.extend(["{dir}/{sample}/prok_bin".format(dir=METABAT_OUTPUT, sample=sample) for sample in sample_df['sample']])

### autometa.smk
if config['autometa']['run_autometa']:
    include: "rules/autometa.smk"

    # input_list.extend(["{dir}/{sample}/intermediates/coverage.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])  # cov_tab
    # input_list.extend(["{dir}/{sample}/intermediates/blastp.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])    # blastp
    # input_list.extend(["{dir}/{sample}/intermediates/taxonomy/taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])   # taxonomy
    
    for kingdom in config['autometa']['binning_target']:
        # input_list.extend(["{dir}/{sample}/intermediates/{kingdom}.markers.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])  # autometa_markers
        # input_list.extend(["{dir}/{sample}/{kingdom}_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])  # binning_output
        # input_list.extend(["{dir}/{sample}/{kingdom}_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # main_output
        if config['autometa']['unclustered_recruitment']:
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_binning.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])    # metabin_stats
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_features.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin_taxonomy
            input_list.extend(["{dir}/{sample}/{kingdom}_recruitment_main.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin
        # input_list.extend(["{dir}/{sample}/{kingdom}_metabin_stats.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']])    # metabin_stats
        input_list.extend(["{dir}/{sample}/{kingdom}_metabin_taxonomy.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin_taxonomy
        # input_list.extend(["{dir}/{sample}/{kingdom}_metabins".format(dir=AUTOMETA_OUTPUT, sample=sample, kingdom=kingdom) for sample in sample_df['sample']]) # metabin
    
    if config['autometa']['align_MAGs']:
        input_list.extend(["{dir}/{sample}/bbmap_refstats.tsv".format(dir=AUTOMETA_OUTPUT, sample=sample) for sample in sample_df['sample']])


############### Rules ######################
localrules: contig_link


wildcard_constraints:
        # ext = "f(ast)?q($|\.gz$)",      # Regex for fastq, fq, fastq.gz and fq.gz as extension
        sample = "[^/]+"                # Regex for all characters except /

rule all:
    input:
        input_list

rule deinterleave:
    input:
        lambda wildcards: os.path.abspath(os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'interleaved'].item()))
    output:
        R1 = "{dir}/{{sample}}_R1.fastq.gz".format(dir=FASTQ),
        R2 = "{dir}/{{sample}}_R2.fastq.gz".format(dir=FASTQ)
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in={input} out1={output.R1} out2={output.R2}
        """

rule fastqc_pre:
    input:
        lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
    output:
        expand("{dir}/pre_trim/{{sample}}_R1_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"]),
        expand("{dir}/pre_trim/{{sample}}_R2_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"])
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
        R1 = lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
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
        rules.trimmomatic.output.R1_paired,
        rules.trimmomatic.output.R2_paired
    output:
        expand("{dir}/post_trim/{{sample}}_R1_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"]),
        expand("{dir}/post_trim/{{sample}}_R2_fastqc.{ext}", dir=FASTQC_OUTPUT, ext=["html", "zip"])
    params:
        dirname = "{dir}/post_trim".format(dir=FASTQC_OUTPUT)
############################################

rule spades:
    input:
        R1 = rules.trimmomatic.output.R1_paired if config['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = rules.trimmomatic.output.R2_paired if config['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
    output:
        "{dir}/{{sample}}/contigs.fasta".format(dir=SPADES_OUTPUT)
    params:
        dirname = directory("{dir}/{{sample}}".format(dir=SPADES_OUTPUT))
    threads: 12
    resources:
        time="14-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (1.5 + attempt * 0.5), 100000), 500000)
        # Set the mem as input_size(mb) * 10 * (2 for first try, 2.5 for second try and 3 for third try) or at least 100G
        # and the maximun usage would not excess 500000 (500G)
    conda:
        "envs/assembler.yaml"
    shell:
        """
        spades.py --meta -o {params.dirname} -1 {input.R1} -2 {input.R2} -t {threads} -m 500
        """

rule megahit:
    input:
        R1 = rules.trimmomatic.output.R1_paired if config['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = rules.trimmomatic.output.R2_paired if config['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
    output:
        "{dir}/{{sample}}/final.contig.fa".format(dir=MEGAHIT_OUTPUT)
    params:
        dirname = directory("{dir}/{{sample}}".format(dir=MEGAHIT_OUTPUT))
    threads: 12
    resources:
        time="14-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (1.5 + attempt * 0.5), 50000), 250000)
    conda:
        "envs/assembler.yaml"
    shell:
        """
        megahit -1 {input.R1} -2 {input.R2} -o {params.dirname} -t {threads}
        """

rule contig_link:
    input:
        rules.spades.output if config['assembler'] == 'spades' else rules.megahit.output
    output:
        "{dir}/{{sample}}_contigs.fasta".format(dir=FILTERED_CONTIGS)
    shell:
        """
        ln -sr {input} {output}
        """

rule filter_contig_length:
    input:
        rules.contig_link.output
    output:
        "{dir}/{{sample}}_filtered.fasta".format(dir=FILTERED_CONTIGS)
    params:
        min_contig_length = config['min_contig_length']
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in={input} out={output} minlength={params.min_contig_length}
        """

rule bbsplit_align_MAGs:
    input:
        genomes = directory("{dir}/{{sample}}/bacteria_metabins".format(dir=AUTOMETA_OUTPUT)),
        R1 = rules.trimmomatic.output.R1_paired if config['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = rules.trimmomatic.output.R2_paired if config['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
    output:
        idx = temp(directory("{dir}/{{sample}}/ref".format(dir=AUTOMETA_OUTPUT))),
        refstats = "{dir}/{{sample}}/bbmap_refstats.tsv".format(dir=AUTOMETA_OUTPUT)
    threads: 8
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        bbsplit.sh -Xmx{resources.mem_mb}m threads={threads} ref={input.genomes} path={output.idx} in={input.R1} in2={input.R2} refstats={output.refstats}
        """



# rule bowtie2_index:
#     input:
#         "{dir}/{{sample}}_filtered.fasta".format(dir=FILTERED_CONTIGS)
#     output:
#         temp(expand("{dir}/{{sample}}.{ext}.bt2", dir=MAPPING_OUTPUT, ext=["1", "2", "3", "4", "rev.1", "rev.2"]))
#     params:
#         idx = "{dir}/{{sample}}".format(dir=MAPPING_OUTPUT)
#     conda:
#         "envs/assembler.yaml"
#     shell:
#         """
#         bowtie2-build {input} {params.idx}
#         """

# rule bowtie2_mapping:
#     input:
#         fastq = lambda wildcards: os.path.join(FASTQ_TRIMMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item()),
#         idx = expand("{dir}/{{sample}}.{ext}.bt2", dir=MAPPING_OUTPUT, ext=["1", "2", "3", "4", "rev.1", "rev.2"])
#     output:
#         bam = "{dir}/{{sample}}.bam".format(dir=MAPPING_OUTPUT),
#         bai = "{dir}/{{sample}}.bam.bai".format(dir=MAPPING_OUTPUT),
#         summary = "{dir}/{{sample}}_align_summary.txt".format(dir=MAPPING_OUTPUT)
#     params:
#         idx = "{dir}/{{sample}}".format(dir=MAPPING_OUTPUT)
#     threads: 8
#     resources:
#         time="1-00:00:00",
#         mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
#     conda:
#         "envs/assembler.yaml"
#     shell:
#         """
#         bowtie2 -p {threads} -x {params.idx} --interleaved {input.fastq} 2> {output.summary} | samtools view -@ {threads} -Sbhu - | samtools sort -@ {threads} -o {output.bam}
#         samtools index {output.bam} {output.bai}
#         """

# rule samtools_idxstats:
#     input:
#         rules.bowtie2_mapping.output.bam
#     output:
#         "{dir}/{{sample}}.stats".format(dir=MAPPING_OUTPUT)
#     shell:
#         """
#         samtools idxstats {input} > {output}
#         """
