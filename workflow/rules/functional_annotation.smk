localrules: link_autometa_prediction

# if config['autometa']['run']:
ruleorder: link_autometa_prediction > prodigal

rule prodigal:
    input:
        "{dir}/{{sample}}_filtered.fasta".format(dir=FILTERED_CONTIGS)
    output:
        nucls = "{dir}/gene_calling/{{sample}}.fna".format(dir=FUNC_ANNO_OUTPUT),
        prots = "{dir}/gene_calling/{{sample}}.faa".format(dir=FUNC_ANNO_OUTPUT)
    params:
        "-p meta -m"
    resources:
        time="7-00:00:00",
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/functional_annotation.yaml"
    shell:
        """
        prodigal -i {input} -a {output.prots} -d {output.nucls} {params}
        """

rule link_autometa_prediction:
    input:
        nucls = "{dir}/{{sample}}/intermediates/predict_orfs.fna".format(dir=AUTOMETA_OUTPUT),
        prots = "{dir}/{{sample}}/intermediates/predict_orfs.faa".format(dir=AUTOMETA_OUTPUT)
    output:
        nucls = "{dir}/gene_calling/{{sample}}.fna".format(dir=FUNC_ANNO_OUTPUT),
        prots = "{dir}/gene_calling/{{sample}}.faa".format(dir=FUNC_ANNO_OUTPUT)
    shell:
        """
        ln -sr {input.nucls} {output.nucls}
        ln -sr {input.prots} {output.prots}
        """

rule cdhit:
    input:
        "{dir}/gene_calling/{{sample}}.faa".format(dir=FUNC_ANNO_OUTPUT)
    output:
        faa = "{dir}/cdhit/{{sample}}.faa".format(dir=FUNC_ANNO_OUTPUT),
        clstr = "{dir}/cdhit/{{sample}}.faa.clstr".format(dir=FUNC_ANNO_OUTPUT)
    threads: 4
    resources:
        time="3-00:00:00",
        mem_mb=lambda w, attempt: 8000 * attempt
    conda:
        "envs/functional_annotation.yaml"
    shell:
        """
        cd-hit -i {input} -o {output} -c 0.95 -M {resources.mem_mb} -T {threads}
        """

rule generate_prediction_bedfile:
    input:
        "{dir}/gene_calling/{{sample}}.fna".format(dir=FUNC_ANNO_OUTPUT)
    output:
        "{dir}/intermediate/{{sample}}.bed".format(dir=FUNC_ANNO_OUTPUT)
    shell:
        """
        awk '$0 ~ /^>/ {{OFS="\t"; sub(/^>/, "", $1); a=$1; sub(/_[0-9]+$/, "", a); print a, $3, $5, $1}}' {input} > {output}
        """

rule bedtools_multicov:
    input:
        bam = "{dir}/{{sample}}.bam".format(dir=MAPPING_OUTPUT),
        bed = rules.generate_prediction_bedfile.output
    output:
        "{dir}/intermediate/{{sample}}_covbed.tsv".format(dir=FUNC_ANNO_OUTPUT)
    conda:
        "envs/functional_annotation.yaml"
    shell:
        """
        multiBamCov -bams {input.bam} -bed {input.bed} > {output}
        """

rule get_alignment_read_counts:
    input:
        clstr = rules.cdhit.output.clstr,
        covbed = rules.bedtools_multicov.output
    output:
        "{dir}/intermediate/{{sample}}_alignment_readcounts.csv".format(dir=FUNC_ANNO_OUTPUT)
    script:
        "scripts/contig_alignment_read_counts.py"

rule dbcanLight_cazyme:
    input:
        rules.cdhit.output.faa if config['functional_annotation']['clustering'] else "{dir}/gene_calling/{{sample}}.faa".format(dir=FUNC_ANNO_OUTPUT)
    output:
        "{dir}/dbcan/{{sample}}/cazymes.tsv".format(dir=FUNC_ANNO_OUTPUT)
    params:
        dir = lambda w, output: os.path.dirname(output[0])
    threads: 8
    resources:
        time="1-00:00:00",
        mem_mb=10000
    conda:
        "envs/functional_annotation.yaml"
    shell:
        """
        dbcanLight -i {input} -o {params.dir} -m cazyme -t {threads}
        """

rule dbcanLight_substrate:
    input:
        rules.cdhit.output.faa if config['functional_annotation']['clustering'] else "{dir}/gene_calling/{{sample}}.faa".format(dir=FUNC_ANNO_OUTPUT)
    output:
        "{dir}/dbcan/{{sample}}/substrates.tsv".format(dir=FUNC_ANNO_OUTPUT)
    params:
        dir = lambda w, output: os.path.dirname(output[0])
    threads: 8
    resources:
        time="1-00:00:00",
        mem_mb=10000
    conda:
        "envs/functional_annotation.yaml"
    shell:
        """
        dbcanLight -i {input} -o {params.dir} -m sub -t {threads}
        """

rule kofamscan:
    input:
        rules.cdhit.output.faa if config['functional_annotation']['clustering'] else "{dir}/gene_calling/{{sample}}.faa".format(dir=FUNC_ANNO_OUTPUT)
    output:
        "{dir}/kofamscan/{{sample}}.txt".format(dir=FUNC_ANNO_OUTPUT)
    params:
        profile = config['functional_annotation']['kofamscan']['profile'],
        ko_list = config['functional_annotation']['kofamscan']['ko_list']
    threads: 4
    resources:
        time="3-00:00:00",
        mem_mb=8000
    conda:
        "envs/functional_annotation.yaml"
    shell:
        """
        exec_annotation -p {params.profile} -k {params.ko_list} --cpu {threads} --tmp-dir {resources.tmpdir} -o {output} -f detail-tsv {input}
        """
