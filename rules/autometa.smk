DATABASES_DIR = config['Path']['Autometa_databases']
ASSEMBLY_OUTPUT = config['Path']['assembly_output']
FILTERED_CONTIGS = config['Path']['filtered_contigs']
AUTOMETA_OUTPUT = config['Path']['autometa_output']
BINNING_INTERMEDIATES = os.path.join(AUTOMETA_OUTPUT, '{sample}', 'intermediates')

### Generate required databases
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

use rule update_marker_database as update_ncbi_database with:
    output:
        nr = ancient("{dir}/ncbi/nr.gz".format(dir=DATABASES_DIR))          # Very time consuming
    params:
        dbdir = directory("{dir}/ncbi".format(dir=DATABASES_DIR)),
        option = "ncbi"
    resources:
        time="7-00:00:00",
        mem_mb=10000

rule diamond:
    input:
        rules.update_ncbi_database.output.nr
    output:
        ancient("{dir}/ncbi/nr.dmnd".format(dir=DATABASES_DIR))             # Very time consuming
    params:
        output = "{}/nr".format(rules.update_ncbi_database.params.dbdir)
    threads: 20
    resources:
        time="7-00:00:00",
        mem_mb=20000
    conda:
        "envs/autometa.yaml"
    shell:
        """
        diamond makedb \
            --in {input} \
            --db {params.output} \
            --threads {threads}
        """
###

rule autometa_length_filter:
    input:
        "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT)
    output:
        fasta = "{dir}/{{sample}}_filtered.fasta".format(dir=FILTERED_CONTIGS),
        gc_content = "{dir}/{{sample}}_filtered_gc_content.tsv".format(dir=FILTERED_CONTIGS)
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
    resources:
        time="7-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
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

rule diamond_blastp:
    input:
        prots_orf = rules.autometa_orf.output.prots,
        db = rules.diamond.output
    output:
        "{dir}/blastp.tsv".format(dir=BINNING_INTERMEDIATES)                # Very time consuming
    threads: 20
    resources:
        time="7-00:00:00",
        mem_mb=20000
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
        "{dir}/lca.tsv".format(dir=BINNING_INTERMEDIATES)                   # time consuming
    conda:
        "envs/autometa.yaml"
    resources:
        time="3-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
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
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    shell:
        """
        autometa-taxonomy-majority-vote \
            --lca {input.lca} \
            --dbdir {input.dbdir} \
            --output {output}
        """

checkpoint autometa_taxonomy:
    input:
        assembly = rules.autometa_length_filter.output.fasta,
        votes = rules.taxonomy_majority_vote.output,
        dbdir = rules.update_ncbi_database.params.dbdir
    output:
        dirname = directory("{dir}/taxonomy".format(dir=BINNING_INTERMEDIATES)),
        taxonomy = "{dir}/taxonomy/taxonomy.tsv".format(dir=BINNING_INTERMEDIATES)
    resources:
        time="7-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-taxonomy \
            --votes {input.votes} \
            --output {output.dirname} \
            --assembly {input.assembly} \
            --ncbi {input.dbdir} \
            --split-rank-and-write superkingdom
        """

def get_superkingdom_fna(wildcards):
    '''
    Retreive the dynamically created files (bacteria.fna and archaea.fna)
    First, retrieve the taxonomy dir from the previous step (checkpoints.autometa_taxonomy)
    Next, directly return the 
    '''
    taxonomy_dir = checkpoints.autometa_taxonomy.get(**wildcards).output[0]
    return "{dir}/{{kingdom}}.fna".format(dir=taxonomy_dir)

rule autometa_kmers:
    input:
        get_superkingdom_fna
    output:
        kmers = "{dir}/{{kingdom}}_kmers.tsv".format(dir=BINNING_INTERMEDIATES),
        kmers_norm = "{dir}/{{kingdom}}_kmers_normalized.tsv".format(dir=BINNING_INTERMEDIATES),
        kmers_embed = "{dir}/{{kingdom}}_kmers_embedded.tsv".format(dir=BINNING_INTERMEDIATES)
    wildcard_constraints:
        kingdom = "bacteria|archaea"
    threads: 20
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-kmers \
            --fasta {input} \
            --kmers {output.kmers} \
            --size 5 \
            --norm-method am_clr \
            --norm-output {output.kmers_norm} \
            --pca-dimensions 50 \
            --embedding-method bhsne \
            --embedding-output {output.kmers_embed} \
            --cpus {threads}
        """

rule autometa_binning:
    input:
        kmers = rules.autometa_kmers.output.kmers_embed,
        cov = rules.autometa_coverage.output,
        gc_content = rules.autometa_length_filter.output.gc_content,
        markers = rules.autometa_markers.output.markers,
        taxonomy = rules.autometa_taxonomy.output.taxonomy
    output:
        binning = "{dir}/{{sample}}/{{kingdom}}_binning.tsv".format(dir=AUTOMETA_OUTPUT),
        main = "{dir}/{{sample}}/{{kingdom}}_main.tsv".format(dir=AUTOMETA_OUTPUT)
    wildcard_constraints:
        kingdom = "bacteria|archaea"
    threads: 10
    resources:
        time="1-00:00:00",
        mem_mb=lambda wildcards, input, attempt: min(max((input.size // 1000000) * 10 * (0.5 + attempt * 0.5), 8000), 250000)
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-binning \
            --kmers {input.kmers} \
            --coverages {input.cov} \
            --gc-content {input.gc_content} \
            --markers {input.markers} \
            --clustering-method dbscan \
            --completeness 20 \
            --purity 90 \
            --cov-stddev-limit 25 \
            --gc-stddev-limit 5 \
            --taxonomy {input.taxonomy} \
            --output-binning {output.binning} \
            --output-main {output.main} \
            --starting-rank superkingdom \
            --rank-filter superkingdom \
            --rank-name-filter {wildcards.kingdom} \
            --cpus {threads}
        """

rule autometa_unclustered_recruitment:
    input:
        kmers = rules.autometa_kmers.output.kmers_norm,
        cov = rules.autometa_coverage.output,
        binning = rules.autometa_binning.output.binning,
        markers = rules.autometa_markers.output.markers,
        taxonomy = rules.autometa_taxonomy.output.taxonomy
    output:
        binning = "{dir}/{{sample}}/{{kingdom}}_recruitment_binning.tsv".format(dir=AUTOMETA_OUTPUT),
        features = "{dir}/{{sample}}/{{kingdom}}_recruitment_features.tsv".format(dir=AUTOMETA_OUTPUT),
        main = "{dir}/{{sample}}/{{kingdom}}_recruitment_main.tsv".format(dir=AUTOMETA_OUTPUT)
    wildcard_constraints:
        kingdom = "bacteria|archaea"
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-unclustered-recruitment \
            --kmers {input.kmers} \
            --coverage {input.cov} \
            --binning {input.binning} \
            --markers {input.markers} \
            --taxonomy {input.taxonomy} \
            --output-binning {output.binning} \
            --output-features {output.features} \
            --output-main {output.main} \
            --classifier random_forest
        """

rule autometa_binning_summary:
    input:
        main = rules.autometa_binning.output.main,
        markers = rules.autometa_markers.output.markers,
        assembly = "{dir}/{{sample}}/scaffolds.fasta".format(dir=ASSEMBLY_OUTPUT),
        dbdir = rules.update_ncbi_database.params.dbdir
    output:
        stats = "{dir}/{{sample}}/{{kingdom}}_metabin_stats.tsv".format(dir=AUTOMETA_OUTPUT),
        taxonomy = "{dir}/{{sample}}/{{kingdom}}_metabin_taxonomy.tsv".format(dir=AUTOMETA_OUTPUT),
        metabins = directory("{dir}/{{sample}}/{{kingdom}}_metabins".format(dir=AUTOMETA_OUTPUT))
    wildcard_constraints:
        kingdom = "bacteria|archaea"
    conda:
        "envs/autometa.yaml"
    shell:
        """
        autometa-binning-summary \
            --binning-main {input.main} \
            --markers {input.markers} \
            --metagenome {input.assembly} \
            --ncbi {input.dbdir} \
            --output-stats {output.stats} \
            --output-taxonomy {output.taxonomy} \
            --output-metabins {output.metabins}
        """