DATA = config['Data']
DATABASES_DIR = config['Autometa_databases']
ASSEMBLY_OUTPUT = os.path.join(DATA, 'spades')
MAPPING_OUTPUT = os.path.join(DATA, 'bowtie2')
BINNING_OUTPUT = os.path.join(DATA, 'autometa')
BINNING_INTERMEDIATES = os.path.join(BINNING_OUTPUT, '{sample}', 'intermediates')

rule autometa_length_filter:
    input:
        "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT)
    output:
        fasta = "{dir}/filtered.fasta".format(dir=BINNING_INTERMEDIATES),
        gc_content = "{dir}/gc_content.tsv".format(dir=BINNING_INTERMEDIATES)
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-length-filter \
            --assembly {input} \
            --cutoff 3000 \
            --output-fasta {output.fasta} \
            --output-gc-content {output.gc_content}
        """

rule autometa_coverage:
    input:
        rules.autometa_length_filter.output.fasta
    output:
        "{dir}/coverage.tsv".format(dir=BINNING_INTERMEDIATES)
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-coverage \
            --assembly {input} \
            --out {output} \
            --from-spades
        """

rule autometa_orf:
    input:
        rules.autometa_length_filter.output.fasta
    output:
        nucls = "{dir}/predict_orfs.fna".format(dir=BINNING_INTERMEDIATES),
        prots = "{dir}/predict_orfs.faa".format(dir=BINNING_INTERMEDIATES)
    threads: 4
    # resources:
    #     mem_mb=240000,
    #     time=20160
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-orfs \
            --assembly {input} \
            --output-nucls {output.nucls} \
            --output-prots {output.prots} \
            --cpus {threads}
        """

rule update_marker_database:
    output:
        hmm_bac = "{dir}/markers/bacteria.single_copy.hmm".format(dir=DATABASES_DIR),
        hmm_arc = "{dir}/markers/archaea.single_copy.hmm".format(dir=DATABASES_DIR)
    params:
        dbdir = directory("{dir}/markers".format(dir=DATABASES_DIR)),
        option = "markers"
    conda:
        "envs/autometa.yaml"
    shell:
        """
        mkdir -p {params.dbdir}

        autometa-config \
            --section databases \
            --option {params.option} \
            --value {params.dbdir}

        autometa-update-databases --update-{params.option}
        """

rule autometa_markers:
    input:
        orfs = rules.autometa_orf.output.prots,
        dbdir = rules.update_marker_database.params.dbdir
    output:
        markers = "{dir}/{{kingdom}}.markers.tsv".format(dir=BINNING_INTERMEDIATES),
        hmmscan = "{dir}/{{kingdom}}.hmmscan.tsv".format(dir=BINNING_INTERMEDIATES)
    wildcard_constraints:
        kingdom = "bacteria|archaea"
    threads: 4
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-markers \
            --orfs {input.orfs} \
            --kingdom {wildcards.kingdom} \
            --dbdir {input.dbdir} \
            --out {output.markers} \
            --hmmscan {output.hmmscan} \
            --parallel \
            --cpus {threads}
        """

use rule update_marker_database as update_ncbi_database with:
    output:
        nr = ancient("{dir}/ncbi/nr.gz".format(dir=DATABASES_DIR))          # Very time consuming
    params:
        dbdir = directory("{dir}/ncbi".format(dir=DATABASES_DIR)),
        option = "ncbi"


rule diamond:
    input:
        rules.update_ncbi_database.output.nr
    output:
        ancient("{dir}/ncbi/nr.dmnd".format(dir=DATABASES_DIR))             # Very time consuming
    params:
        output = "{}/nr".format(rules.update_ncbi_database.params.dbdir)
    threads: 20
    conda:
        "envs/autometa.yaml"
    shell:
        """
        diamond makedb \
            --in {input} \
            --db {params.output} \
            --threads {threads}
        """

rule diamond_blastp:
    input:
        prots_orf = rules.autometa_orf.output.prots,
        db = rules.diamond.output
    output:
        "{dir}/blastp.tsv".format(dir=BINNING_INTERMEDIATES)
    threads: 4
    conda:
        "envs/autometa.yaml"
    shell:
        """
        diamond blastp \
            --query {input.prots_orf} \
            --db {input.db} \
            --evalue 1e-5 \
            --max-target-seqs 200 \
            --threads {threads} \
            --outfmt 6 \
            --out {output}
        """

rule autometa_taxonomy_lca:
    input:
        blastp = rules.diamond_blastp.output,
        dbdir = rules.update_ncbi_database.params.dbdir
    output:
        "{dir}/lca.tsv".format(dir=BINNING_INTERMEDIATES)
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-taxonomy-lca \
            --blast {input.blastp} \
            --dbdir {input.dbdir} \
            --lca-output {output}
        """

rule taxonomy_majority_vote:
    input:
        lca = rules.autometa_taxonomy_lca.output,
        dbdir = rules.update_ncbi_database.params.dbdir
    output:
        "{dir}/votes.tsv".format(dir=BINNING_INTERMEDIATES)
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-taxonomy-majority-vote \
            --lca {input.lca} \
            --dbdir {input.dbdir} \
            --output {output}
        """

rule autometa_taxonomy:
    input:
        assembly = rules.autometa_length_filter.output.fasta,
        votes = rules.taxonomy_majority_vote.output,
        dbdir = rules.update_ncbi_database.params.dbdir
    output:
        directory("{dir}/taxonomy".format(dir=BINNING_INTERMEDIATES))
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-taxonomy \
            --votes {input.votes} \
            --output {output} \
            --assembly {input.assembly} \
            --ncbi {input.dbdir} \
            --split-rank-and-write superkingdom
        """
