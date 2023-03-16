GTDBTK_DB = config['gtdbtk']['db']
MASH_DB = config['gtdbtk']['mash_db']

rule gtdbtk_classify_wf:
    input:
        "{dir}/{{sample}}".format(dir=METABAT_OUTPUT),
    output:
        directory("{dir}/{{sample}}".format(dir=GTDBTK_OUTPUT))
    threads: 2
    params:
        mash = MASH_DB
    resources:
        time="7-00:00:00",
        mem_mb=80000
    conda:
        "envs/gtdbtk.yaml"
    shell:
        """
        gtdbtk classify_wf --genome_dir {input} --out_dir {output} -x fa --cpus {threads} --keep_intermediates --mash_db {params.mash}
        """
