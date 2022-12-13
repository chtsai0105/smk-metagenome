localrules: contig_link

if config['assembly']['assembler'] == 'spades':
    ruleorder: spades > megahit
else:
    ruleorder: megahit > spades


rule spades:
    input:
        R1 = rules.trimmomatic.output.R1_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = rules.trimmomatic.output.R2_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
    output:
        "{dir}/{{sample}}/contigs.fasta".format(dir=ASSEMBLY_OUTPUT)
    params:
        dirname = directory("{dir}/{{sample}}".format(dir=ASSEMBLY_OUTPUT))
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
        R1 = rules.trimmomatic.output.R1_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = rules.trimmomatic.output.R2_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
    output:
        "{dir}/{{sample}}/final.contig.fa".format(dir=ASSEMBLY_OUTPUT)
    params:
        dirname = directory("{dir}/{{sample}}".format(dir=ASSEMBLY_OUTPUT))
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
        rules.spades.output if config['assembly']['assembler'] == 'spades' else rules.megahit.output
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
        min_contig_length = config['assembly']['min_contig_length']
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in={input} out={output} minlength={params.min_contig_length}
        """
