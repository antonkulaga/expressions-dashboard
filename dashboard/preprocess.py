import genotations.genomes
from pathlib import Path
from pycomfort.files import *
import pyarrow
import pandas as pd
from functional import seq
from typing import *
import functools


import polars as pl

pl.col("some_string_column").str.split(".")
pl.col("some_array_column").arr.sum()

def read_quant(path: Path, transcripts: bool) -> pl.DataFrame:
    alias = "transcript" if transcripts else "gene"
    gene = pl.col("Name").str.split(".").apply(lambda s: s[0]).alias(alias)
    dtypes={"TPM": pl.datatypes.Float64, "EffectiveLength": pl.datatypes.Float64, "NumReads": pl.datatypes.Float64}
    return pl.read_csv(path, sep="\t", dtypes = dtypes).with_column(gene).select([gene, pl.col("TPM"), pl.col("EffectiveLength"), pl.col("NumReads")])

def quants_from_bioproject(project: Path, name_part: str = "quant.sf", dir_part: str = "quant_") -> OrderedDict[str, pl.DataFrame]:
    """
    reads Salmon quant.sf files from the folder and writes them into Ordered dictionary of polars dataframes
    :param project:
    :param name_part:
    :param dir_part:
    :return:
    """
    transcripts = "gene" not in dir_part and "gene" not in name_part
    return OrderedDict([
        (
            q.name.replace(dir_part, ""),
            read_quant(files(q).filter(lambda f: name_part in f.name).first(), transcripts)
        )
        for q in dirs(project).flat_map(lambda run: dirs(run).filter(lambda f: dir_part in f.name))
    ])


def with_mouse_transcript_info(expressions: pl.DataFrame):
    return genotations.genomes.mouse.annotations.transcript_gene_names_df.join(expressions, on="transcript")

def with_mouse_gene_info(expressions: pl.DataFrame):
    return genotations.genomes.mouse.annotations.gene_names_df.join(expressions, on="gene")


def transcripts_from_bioproject(project: Path) -> OrderedDict[str, pl.DataFrame]:
    return quants_from_bioproject(project, "quant.sf")


def genes_from_bioproject(project: Path) -> OrderedDict[str, pl.DataFrame]:
    return quants_from_bioproject(project, "quant.genes.sf")


def expressions_from_bioproject(project: Path, transcripts: bool = True) -> pl.DataFrame:
    """
    Merges expression from multiple Salmon quant files to one expressions/samples matrix
    :param project:
    :param transcripts:
    :return:
    """
    expressions: OrderedDict = transcripts_from_bioproject(project) if transcripts else genes_from_bioproject(project)
    name = "transcript" if transcripts else "gene"
    frames = [v.select([pl.col(name), pl.col("TPM").alias(k)]) for k, v in expressions.items()]
    return functools.reduce(lambda a, b: a.join(b, on=name), frames)