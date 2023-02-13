DB_PATH = config['kraken2']['db_dir']
DB = config['kraken2']['db_name']


rule kraken2_database:
    output:
        "{dir}/{db}".format(dir=DB_PATH, db=DB)        # Extremly time consuming
    threads: 12
    params:
        db_path = DB_PATH,
        db = DB
    resources:
        time="7-00:00:00",
        mem_mb=120000
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        KRAKEN2_DB_PATH={params.db_path}
        kraken2-build --download-taxonomy --db {params.db} --threads {threads}
        kraken2-build --download-library fungi --db {params.db} --threads {threads}
        kraken2-build --download-library bacteria --db {params.db} --threads {threads}
        kraken2-build --build --db {params.db}
        kraken2-build --clean --db {params.db}
        """

rule kraken2_profiling:
    input:
        R1 = lambda w: data.unified_samples(w.sample, FASTQ_TRIMMED, read='R1', ext='.fastq.gz') if config['trimming']['run']
            else data.unified_samples(w.sample, FASTQ, read='R1', column='renamed_fastq'),
        R2 = lambda w: data.unified_samples(w.sample, FASTQ_TRIMMED, read='R2', ext='.fastq.gz') if config['trimming']['run']
            else data.unified_samples(w.sample, FASTQ, read='R2', column='renamed_fastq'),
        database = rules.kraken2_database.output
    output:
        "{dir}/{{sample}}_kraken2.tsv".format(dir=KRAKEN2_OUTPUT)
    threads: 8
    params:
        db_path = DB_PATH,
        db = DB
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        KRAKEN2_DB_PATH={params.db_path}
        kraken2 --db {params.db_path}/{params.db} --threads {threads} --paired --gzip-compressed --confidence 0.1 {input.R1} {input.R2} --output {output}
        """

rule taxonkit_reformat:
    input:
        rules.kraken2_profiling.output
    output:
        taxid = temp("{dir}/{{sample}}_taxid.tsv".format(dir=KRAKEN2_OUTPUT)),
        taxon_csv = "{dir}/{{sample}}_taxon.csv".format(dir=KRAKEN2_OUTPUT)
    conda:
        "envs/kraken2.yaml"
    shell:
        """
        cut -f3 {input} | sort | uniq -c | sort -nr | awk '{{print $2"\t"$1}}' > {output.taxid}
        taxonkit reformat -I 1 -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" {output.taxid} > {output.taxon_csv}
        sed -i '1 i\ttaxid\treads\tsuperkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies' {output.taxon_csv}
        sed -i 's/\t/,/g' {output.taxon_csv}
        """