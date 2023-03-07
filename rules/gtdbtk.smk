rule gtdbtk_identify:
    input:
        "{dir}/{{sample}}/bin".format(dir=METABAT_OUTPUT)
    output:
        directory("{dir}/{{sample}}/identify".format(dir=GTDBTK_OUTPUT))
    threads: 2
    resources:
        time="3-00:00:00",
        mem_mb=20000
    conda:
        "envs/gtdbtk.yaml"
    shell:
        """
        gtdbtk identify --genome_dir {input} --out_dir {output} --extension fasta --cpus {threads}
        """

rule gtdbtk_align:
    input:
        rules.gtdbtk_identify.output
    output:
        directory("{dir}/{{sample}}/align".format(dir=GTDBTK_OUTPUT))
    threads: 2
    resources:
        time="3-00:00:00",
        mem_mb=20000
    conda:
        "envs/gtdbtk.yaml"
    shell:
        """
        gtdbtk align --identify_dir {input} --out_dir {output} --cpus {threads}
        """

rule gtdbtk_classify:
    input:
        fasta = "{dir}/{{sample}}/bin".format(dir=METABAT_OUTPUT),
        align = rules.gtdbtk_align.output
    output:
        directory("{dir}/{{sample}}/classify".format(dir=GTDBTK_OUTPUT))
    threads: 2
    resources:
        time="7-00:00:00",
        mem_mb=100000
    conda:
        "envs/gtdbtk.yaml"
    shell:
        """
        gtdbtk classify --genome_dir {input.fasta} --align_dir {input.align} --out_dir {output} -x fasta --cpus {threads}
        """
