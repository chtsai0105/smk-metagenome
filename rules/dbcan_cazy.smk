rule run_dbcan:
    input:
        "{dir}/{{sample}}/intermediates/predict_orfs.faa".format(dir=AUTOMETA_OUTPUT)
    output:
        dirname = directory("{dir}/{{sample}}".format(dir=DBCAN_OUTPUT)),
        hmmer = "{dir}/{{sample}}/hmmer.out".format(dir=DBCAN_OUTPUT),
        overview = "{dir}/{{sample}}/overview.txt".format(dir=DBCAN_OUTPUT)
    params:
        db = config['dbcan_cazy']['db']
    threads: 8
    resources:
        time="7-00:00:00",
        mem_mb=80000
    conda:
        "envs/dbcan.yaml"
    shell:
        """
        run_dbcan --db_dir {params.db} --out_dir {output.dirname} --tools hmmer \
        --hmm_cpu {threads} --dia_cpu {threads} --dbcan_thread {threads} --tf_cpu {threads} --stp_cpu {threads} \
        {input} protein
        """
