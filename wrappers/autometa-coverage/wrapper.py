"""Snakemake wrapper for Autometa subprocess coverage."""

__author__ = "Cheng-Hung Tsai @ UCR"
__copyright__ = "Copyright 2023, Cheng-Hung Tsai"
__email__ = "chenghung.tsai@email.ucr.edu"
__license__ = "MIT"

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

if hasattr(snakemake.params, "spades") and snakemake.params.get("spades"):
    shell(
        "autometa-coverage "
        " --assembly {snakemake.input.fasta} "
        " --out {snakemake.output} "
        " --from-spades"
    )

else:
    if not hasattr(snakemake.input, "bam") or not snakemake.input.get("bam"):
        raise KeyError("When assembled by other assembler than spades, alignment bam file should be specified.")

    shell(
        "autometa-coverage "
        " --assembly {snakemake.input.fasta} "
        " --bam {snakemake.input.bam} "
        " --out {snakemake.output} "
        " --cpus {snakemake.threads}"
    )
