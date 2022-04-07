DATA = config['Data']
ASSEMBLY_OUTPUT = os.path.join(DATA, 'spades')
MAPPING_OUTPUT = os.path.join(DATA, 'bowtie2')
EUKREP_OUTPUT = os.path.join(DATA, 'EukRep')
EUK_BINNING_OUTPUT = os.path.join(DATA, 'metabat_euk')

rule eukrep:
    input:
        "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT)
    output:
        "{dir}/{{sample}}_euk.fasta".format(dir=EUKREP_OUTPUT)
    conda:
        "envs/euk_detection.yaml"
    shell:
        """
        EukRep -i {input} -o {output} --min 1000
        """

rule metabat_contig_depths:
    input:
        bam = "{dir}/{{sample}}.bam".format(dir=MAPPING_OUTPUT),
        bai = "{dir}/{{sample}}.bam.bai".format(dir=MAPPING_OUTPUT)
    output:
        "{dir}/{{sample}}_cov_metabat.tsv".format(dir=MAPPING_OUTPUT),
    conda:
        "envs/euk_detection.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam}
        """

rule metabat_binning_euk:
    input:
        assembly = rules.eukrep.output,
        cov_tab = rules.metabat_contig_depths.output
    output:
        "{dir}/{{sample}}/euk_bin".format(dir=EUK_BINNING_OUTPUT)
    threads: 4
    conda:
         "envs/euk_detection.yaml"
    shell:
        """
        metabat2 --numThreads {threads} --saveCls -i {input.assembly} -a {input.cov_tab} -o {output}
        """