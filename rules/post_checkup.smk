localrules: add_unmapped_reads

rule bbsplit_align_MAGs:
    input:
        genomes = directory("{dir}/{{sample}}/bacteria_metabins".format(dir=AUTOMETA_OUTPUT)),
        R1 = rules.trimmomatic.output.R1_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = rules.trimmomatic.output.R2_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item())
    output:
        idx = temp(directory("{dir}/{{sample}}_ref".format(dir=MAPPING_OUTPUT))),
        refstats = temp("{dir}/{{sample}}_temp_refstats.tsv".format(dir=MAPPING_OUTPUT)),
        log = "{dir}/{{sample}}.log".format(dir=MAPPING_OUTPUT)
    threads: 8
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        bbsplit.sh -Xmx{resources.mem_mb}m threads={threads} ref={input.genomes} path={output.idx} in={input.R1} in2={input.R2} refstats={output.refstats} 2> {output.log}
        """

rule add_unmapped_reads:
    input:
        refstats = rules.bbsplit_align_MAGs.output.refstats,
        log = rules.bbsplit_align_MAGs.output.log
    output:
        "{dir}/{{sample}}_refstats.tsv".format(dir=MAPPING_OUTPUT)
    shell:
        """
        total_reads=`awk '/Reads\ Used/{{gsub(/[\ ]/,""); print$2}}' {input.log}`
        unmapped_reads=$(expr $total_reads - `awk '{{read+=$8}} END {{print read}}' {input.refstats}`)
        echo -e "unmapped\t*\t*\t*\t*\t*\t*\t$unmapped_reads\t*" >> {input.refstats}
        cp {input.refstats} {output}
        """

rule bowtie2_index:
    input:
        "{dir}/{{sample}}_filtered.fasta".format(dir=FILTERED_CONTIGS)
    output:
        temp(expand("{dir}/{{sample}}.{ext}.bt2", dir=MAPPING_OUTPUT, ext=["1", "2", "3", "4", "rev.1", "rev.2"]))
    params:
        idx = "{dir}/{{sample}}".format(dir=MAPPING_OUTPUT)
    conda:
        "envs/assembler.yaml"
    shell:
        """
        bowtie2-build {input} {params.idx}
        """

rule bowtie2_mapping:
    input:
        R1 = rules.trimmomatic.output.R1_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R1'].item()),
        R2 = rules.trimmomatic.output.R2_paired if config['trimming']['run_trimmomatic'] else lambda wildcards: os.path.join(FASTQ, sample_df.loc[sample_df['sample'] == wildcards.sample, 'R2'].item()),
        idx = expand("{dir}/{{sample}}.{ext}.bt2", dir=MAPPING_OUTPUT, ext=["1", "2", "3", "4", "rev.1", "rev.2"])
    output:
        bam = "{dir}/{{sample}}.bam".format(dir=MAPPING_OUTPUT),
        bai = "{dir}/{{sample}}.bam.bai".format(dir=MAPPING_OUTPUT),
        summary = "{dir}/{{sample}}_align_summary.txt".format(dir=MAPPING_OUTPUT)
    params:
        idx = "{dir}/{{sample}}".format(dir=MAPPING_OUTPUT)
    threads: 8
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/assembler.yaml"
    shell:
        """
        bowtie2 -p {threads} -x {params.idx} -1 {input.R1} -2 {input.R2} 2> {output.summary} | samtools view -@ {threads} -Sbhu - | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam} {output.bai}
        """

rule samtools_idxstats:
    input:
        rules.bowtie2_mapping.output.bam
    output:
        "{dir}/{{sample}}.stats".format(dir=MAPPING_OUTPUT)
    shell:
        """
        samtools idxstats {input} > {output}
        """