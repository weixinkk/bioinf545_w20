"""Snakemake workflow for b545 project
"""
from pathlib import Path
import pandas as pd
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

configfile: "config.yaml"
conda: "environment.yml"
int_dir = Path(config["intermediate_dir"])
out_dir = Path(config["output_dir"])

manifest = pd.read_csv(config["manifest"],sep="\t")

FTP = FTPRemoteProvider()

rule all:
    input: expand(str(int_dir / "quant" / "{SRR_ID}" / "abundance.tsv"),SRR_ID=manifest["SRR"])

rule fastq_dump:
    output: 
        R1=str(int_dir / "fastqs" / "{SRR_ID}" / "{SRR_ID}_1.fastq.gz"),
        R2=str(int_dir / "fastqs" / "{SRR_ID}" / "{SRR_ID}_2.fastq.gz")
    shell:
        """
        mkdir -p $(dirname {output.R1})
        fastq-dump -O $(dirname {output.R1}) --split-files --gzip {wildcards.SRR_ID}
        """

rule kallisto_make_index:
    input: FTP.remote(config["reference"])
    output: str(int_dir/"reference_index.idx")
    shell:
        """
        kallisto index -i {output} {input}
        """
    
rule kallisto_quant:
    input: 
        R1=str(int_dir / "fastqs" / "{SRR_ID}" / "{SRR_ID}_1.fastq.gz"),
        R2=str(int_dir / "fastqs" / "{SRR_ID}" / "{SRR_ID}_2.fastq.gz"),
        index=str(int_dir/"reference_index.idx")
    output:
        h5=str(int_dir / "quant" / "{SRR_ID}" / "abundance.h5"),
        abundance=str(int_dir / "quant" / "{SRR_ID}" / "abundance.tsv"),
        run_info=str(int_dir / "quant" / "{SRR_ID}" / "run_info.json")
    shell:
        """
        mkdir -p $(dirname {output.h5})
        kallisto quant -i {input.index} -o $(dirname {output.h5}) {input.R1} {input.R2}
        """
