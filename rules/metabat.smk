rule metabat_contig_depths:
    input:
        bam = "{dir}/{{sample}}.bam".format(dir=MAPPING_OUTPUT),
        bai = "{dir}/{{sample}}.bam.bai".format(dir=MAPPING_OUTPUT)
    output:
        "{dir}/{{sample}}_cov_metabat.tsv".format(dir=METABAT_OUTPUT),
    conda:
        "envs/metabat.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam}
        """

rule metabat_binning:
    input:
        fasta = "{dir}/{{sample}}_filtered.fasta".format(dir=FILTERED_CONTIGS),
        cov = rules.metabat_contig_depths.output
    output:
        "{dir}/{{sample}}/bin".format(dir=METABAT_OUTPUT)
    threads: 4
    resources:
        time="3-00:00:00",
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/metabat.yaml"
    shell:
        """
        metabat2 --saveCls --numThreads {threads} -i {input.fasta} -a {input.cov} -o {output}
        """