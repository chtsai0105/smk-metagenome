rule mmseq_easy_taxonomy:
    input:
        R1 = "{dir}/preprocess_done/{{sample}}_R1.fastq.gz".format(dir=FASTQ),
        R2 = "{dir}/preprocess_done/{{sample}}_R2.fastq.gz".format(dir=FASTQ),
        db = config['mmseqs2']['db']
    output:
        directory("{dir}/{{sample}}".format(dir=config['mmseqs2']['output']))
    threads: 8
    resources:
        time="3-00:00:00",
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * (1 + attempt * 0.5), 40000), 250000)
    conda:
        "envs/mmseqs2.yaml"
    shell:
        """
        mmseqs easy-taxonomy {input.R1} {input.R2} {input.db} {output} {resources.tmpdir} \
        --threads {threads} --tax-lineage 1 --lca-ranks kingdom,phylum,family,genus --db-load-mode 2 \
        --split-memory-limit {resources.mem_mb}
        """
