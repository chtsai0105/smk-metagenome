localrules: contig_link

if config['assembly']['assembler'] == 'spades':
    ruleorder: spades > megahit
else:
    ruleorder: megahit > spades

rule spades:
    input:
        reads = expand("{dir}/preprocess_done/{{sample}}_{R}.fastq.gz", dir=FASTQ, R=["R1", "R2"])
    output:
        contigs = "{dir}/{{sample}}/contigs.fasta".format(dir=ASSEMBLY_OUTPUT),
        scaffolds = "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT)
    log:
        "logs/spades_{sample}.log"
    params:
        extra="--only-assembler"
    threads: 12
    resources:
        partition=config['partition']['highmem'],
        time="30-00:00:00",
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * 10 * (2 + attempt), 100000), 600000)
        # Set the mem as input_size(mb) * 10 * (3 for first try, 4 for second try and 5 for third try) or at least 100G
        # and the maximun usage would not excess 600000 (600G)
    conda:
        "envs/assembler.yaml"
    wrapper:
        "file:wrappers/metaspades"

rule megahit:
    input:
        R1 = "{dir}/preprocess_done/{{sample}}_R1.fastq.gz".format(dir=FASTQ),
        R2 = "{dir}/preprocess_done/{{sample}}_R2.fastq.gz".format(dir=FASTQ)
    output:
        "{dir}/{{sample}}/final.contigs.fa".format(dir=ASSEMBLY_OUTPUT)
    params:
        dir = lambda w, output: os.path.dirname(output[0])
    threads: 12
    resources:
        time="14-00:00:00",
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * 10 * (1.5 + attempt * 0.5), 50000), 250000)
    conda:
        "envs/assembler.yaml"
    shell:
        """
        rm {params.dir} -rf
        megahit -1 {input.R1} -2 {input.R2} -o {params.dir} -t {threads}
        """

rule contig_link:
    input:
        rules.spades.output.contigs if config['assembly']['assembler'] == 'spades' else rules.megahit.output
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
