FILTERED_CONTIGS = config['Path']['filtered_contigs']
MAPPING_OUTPUT = config['Path']['mapping_output']
EUKREP_OUTPUT = config['Path']['eukrep_output']
METABAT_OUTPUT = config['Path']['metabat_output']

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

rule metabat_binning:
    input:
        assembly = "{dir}/{{sample}}_filtered.fasta".format(dir=MAPPING_OUTPUT),
        cov_tab = rules.metabat_contig_depths.output
    output:
        "{dir}/{{sample}}/bin".format(dir=METABAT_OUTPUT)
    threads: 5
    resources:
        mem_mb=lambda wildcards, input, attempt: max((input.size // 1000000) * 500, 20000) * attempt
    conda:
        "envs/euk_detection.yaml"
    shell:
        """
        metabat2 --saveCls -m 3000 --numThreads {threads} -i {input.assembly} -a {input.cov_tab} -v -o {output}
        """

rule eukrep:
    input:
        "{dir}/{{sample}}_filtered.fasta".format(dir=MAPPING_OUTPUT)
    output:
        euk = "{dir}/{{sample}}_euk.fasta".format(dir=EUKREP_OUTPUT),
        prok = "{dir}/{{sample}}_prok.fasta".format(dir=EUKREP_OUTPUT)
    resources:
        mem_mb=lambda wildcards, input, attempt: max((input.size // 1000000) * 500, 20000) * attempt
    conda:
        "envs/euk_detection.yaml"
    shell:
        """
        EukRep -i {input} -o {output.euk} --prokarya {output.prok} --min 3000
        """

use rule metabat_binning as metabat_binning_euk with:
    input:
        assembly = rules.eukrep.output.euk,
        cov_tab = rules.metabat_contig_depths.output
    output:
        "{dir}/{{sample}}/euk_bin".format(dir=METABAT_OUTPUT)

use rule metabat_binning as metabat_binning_prok with:
    input:
        assembly = rules.eukrep.output.prok,
        cov_tab = rules.metabat_contig_depths.output
    output:
        "{dir}/{{sample}}/prok_bin".format(dir=METABAT_OUTPUT)
