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

DATA = config['Data']
FASTQ = os.path.join(DATA, 'fastq')
FASTQ_RENAMED = os.path.join(DATA, 'fastq_renamed')
FASTQ_TRIMMED = os.path.join(DATA, 'fastq_trimmed')
FASTQC_OUTPUT = os.path.join(DATA, 'fastqc')

ASSEMBLY_OUTPUT = os.path.join(DATA, 'spades')
MAPPING_OUTPUT = os.path.join(DATA, 'bowtie2')
BINNING_OUTPUT = os.path.join(DATA, 'autometa')

AUTOMETA_DATABASES = "/srv/projects/db/autometa/1.0.2"

localrules: rename_input
wildcard_constraints:
        ext = "f(ast)?q($|\.gz$)",      # Regex for fastq, fq, fastq.gz and fq.gz as extension
        sample = "[^/]+"                # Regex for all characters except /

rule all:
    input:
        fastq_renamed =     expand("{dir}/{fastq}", dir=FASTQ_RENAMED, fastq=sample_df['fastq_renamed']),
        fastqc_pre_trim =   expand("{dir}/pre_trim/{sample}_fastqc.{ext}", dir=FASTQC_OUTPUT, sample=sample_df['sample'], ext=["html", "zip"]),
        fastq_trimmed =     expand("{dir}/{fastq}", dir=FASTQ_TRIMMED, fastq=sample_df['fastq_renamed']),
        fastqc_post_trim =  expand("{dir}/post_trim/{sample}_fastqc.{ext}", dir=FASTQC_OUTPUT, sample=sample_df['sample'], ext=["html", "zip"]),
        assembly =          expand("{dir}/{sample}/scaffolds.fasta", dir=ASSEMBLY_OUTPUT, sample=sample_df['sample']),
        cov_tab =           expand("{dir}/{sample}/coverage.tab", dir=MAPPING_OUTPUT, sample=sample_df['sample']),
        # binning =           expand("{dir}/{sample}/{kingdom}_run", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['Bacteria', 'Archaea']),
        cluster_output =    expand("{dir}/{sample}/{kingdom}_run/cluster_process_output", dir=BINNING_OUTPUT, sample=sample_df['sample'], kingdom=['Bacteria', 'Archaea'])

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
    threads: 2
    conda:
        "envs/preprocess.yaml"
    params:
        dirname = "{dir}/pre_trim".format(dir=FASTQC_OUTPUT)
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
        lambda wildcards: os.path.join(FASTQ_TRIMMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
    output:
        dirname = directory("{dir}/{{sample}}".format(dir=ASSEMBLY_OUTPUT)),
        assembly = "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT)
    threads: 12
    resources:
        mem_mb=240
    conda:
        "envs/assembler.yaml"
    shell:
        """
        spades.py --meta -o {output.dirname} --12 {input} -t {threads}
        """

rule bowtie2_mapping:
    input:
        assembly = rules.metaspades.output.assembly,
        fastq = lambda wildcards: os.path.join(FASTQ_TRIMMED, sample_df.loc[sample_df['sample'] == wildcards.sample, 'fastq_renamed'].item())
    output:
        bam = temp("{dir}/{{sample}}/mapped.bam".format(dir=MAPPING_OUTPUT)),
        bai = temp("{dir}/{{sample}}/mapped.bam.bai".format(dir=MAPPING_OUTPUT))
    threads: 4
    envmodules:
        "autometa/1.0.2"
    params:
        idx = "{dir}/{{sample}}/idx".format(dir=MAPPING_OUTPUT)
    shell:
        """
        bowtie2-build {input.assembly} {params.idx}
        bowtie2 -x {params.idx} --interleaved {input.fastq} | samtools -Sbh -o {output}
        samtools index {output}
        rm {params.idx}.*.bt2
        """

rule generate_coverage_table:
    input:
        bam = rules.bowtie2_mapping.output.bam,
        bai = rules.bowtie2_mapping.output.bai
    output:
        cov_bed = temp("{dir}/{{sample}}/genome_cov.bed".format(dir=MAPPING_OUTPUT)),
        cov_tab = "{dir}/{{sample}}/coverage.tab".format(dir=MAPPING_OUTPUT)
    threads: 4
    shell:
        """
        genomeCoverageBed -ibam {input.bam} > {output.cov_bed}
        contig_coverage_from_bedtools.pl {output.cov_bed} > {output.cov_tab}
        """

rule autometa_split_to_kingdom_bins:
    input:
        assembly = rules.metaspades.output.assembly,
        cov_tab = rules.generate_coverage_table.output.cov_tab
    output:
        dirpath = directory("{dir}/{{sample}}".format(dir=BINNING_OUTPUT)),
        fasta = expand("{dir}/{{sample}}/{kingdom}.fasta", dir=BINNING_OUTPUT, kingdom=["Bacteria", "Archaea"]),
        taxon_tab = "{dir}/{{sample}}/taxonomy.tab".format(dir=BINNING_OUTPUT)
    threads: 4
    shell:
        """
        make_taxonomy_table.py \
        --assembly {input.assembly} \
        --processors {threads} \
        --length_cutoff 3000 \
        --output_dir {output.dirpath} \
        --cov_table {input.cov_tab}
        """

rule autometa_binning_by_kingdom:
    input:
        fasta = "{dir}/{{sample}}/{{kingdom}}.fasta".format(dir=BINNING_OUTPUT),
        taxon_tab = rules.autometa_split_to_kingdom_bins.output.taxon_tab,
        cov_tab = rules.generate_coverage_table.output.cov_tab
    output:
        dirpath = directory("{dir}/{{sample}}/{{kingdom}}_run".format(dir=BINNING_OUTPUT)),
        bin_tab = "{dir}/{{sample}}/{{kingdom}}_run/ML_recruitment_output.tab".format(dir=BINNING_OUTPUT)
    params:
        kingdom = lambda wildcards: wildcards.kingdom.lower()
    wildcard_constraints:
        kingdom = "Bacteria|Archaea"
    shell:
        """
        run_autometa.py \
        --kingdom {params.kingdom} \
        --assembly {input.fasta} \
        --processors {threads} \
        --length_cutoff 1500 \
        --taxonomy_table {input.taxon_tab} \
        --output_dir {output.dirpath} \
        --cov_table {input.cov_tab}
        --ML_recruitment \
        """

rule autometa_cluster_analysis:
    input:
        fasta = "{dir}/{{sample}}/{{kingdom}}.fasta".format(dir=BINNING_OUTPUT),
        bin_tab = rules.autometa_binning_by_kingdom.output.bin_tab
    output:
        directory("{dir}/{{sample}}/{{kingdom}}_run/cluster_process_output".format(dir=BINNING_OUTPUT))
    params:
        databases = AUTOMETA_DATABASES
    shell:
        """
        cluster_process.py \
        --bin_table {input.bin_tab} \
        --column ML_expanded_clustering \
        --fasta {input.fasta} \
        --do_taxonomy \
        --db_dir {params.databases} \
        --output_dir {output}
        """
