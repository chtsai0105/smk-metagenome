wildcard_constraints:
        sample = "[^/]+",                # Regex for all characters except /
        R = "R(1|2)",
        ext = "f(ast)?q($|\.gz$)"

rule convert_paired_to_interleaved:
    input:
        R1 = lambda w: data.filepath_generator(w.sample, FASTQ, column="R1"),
        R2 = lambda w: data.filepath_generator(w.sample, FASTQ, column="R2")
    output:
        "{dir}/{{sample}}_interleaved.fastq.gz".format(dir=FASTQ)
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in1={input.R1} in2={input.R2} out={output}
        """

rule fastqc_pre:
    input:
        lambda w: data.filepath_generator(w.sample, FASTQ)
    output:
        html = "{dir}/pre_trim/{{sample}}_fastqc.html".format(dir=FASTQC_OUTPUT),
        zip = "{dir}/pre_trim/{{sample}}_fastqc.zip".format(dir=FASTQC_OUTPUT)
    params:
        dir = lambda w, output: os.path.dirname(output[0]),
        basename = lambda w, output: data.filepath_generator(w.sample, column="base")
    threads: 4
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        fastqc -t {threads} -o {params.dir} {input}
        if [ "{params.dir}/{params.basename}_fastqc.zip" != "{output.zip}" ]; then
            mv {params.dir}/{params.basename}_fastqc.html {output.html}
            mv {params.dir}/{params.basename}_fastqc.zip {output.zip}
        fi
        """

rule fastp:
    input:
        lambda w: data.filepath_generator(w.sample, FASTQ)
    output:
        "{dir}/trimmed/{{sample}}.fastq.gz".format(dir=FASTQ)
    threads: 4
    resources:
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * 4 * (2 + attempt), 4000), 16000)
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        fastp -i {input} -w {threads} --interleaved_in --stdout --length_required 75 | gzip -1 > {output}
        """

use rule fastqc_pre as fastqc_post with:
    input:
        rules.fastp.output
    output:
        html = "{dir}/post_trim/{{sample}}_fastqc.html".format(dir=FASTQC_OUTPUT),
        zip = "{dir}/post_trim/{{sample}}_fastqc.zip".format(dir=FASTQC_OUTPUT)

rule error_correction:
    input:
        rules.fastp.output if config['trimming']['run'] else lambda w: data.filepath_generator(w.sample, FASTQ)
    output:
        "{dir}/corrected/{{sample}}.fastq.gz".format(dir=FASTQ)
    threads: 8
    resources:
        mem_mb=lambda w, input, attempt: min(max((input.size // 1000000) * 10 * (2 + attempt), 20000), 100000)
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        bfc -1 -k 21 -t {threads} {input} | gzip -1 > {output}
        """

rule deinterleave:
    input:
        rules.error_correction.output
    output:
        R1 = "{dir}/preprocess_done/{{sample}}_R1.fastq.gz".format(dir=FASTQ),
        R2 = "{dir}/preprocess_done/{{sample}}_R2.fastq.gz".format(dir=FASTQ)
    params:
        min_read_length = 75
    conda:
        "envs/preprocess.yaml"
    shell:
        """
        reformat.sh in={input} out1={output.R1} out2={output.R2} minlength={params.min_read_length}
        """
