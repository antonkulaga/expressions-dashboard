#!/usr/bin/env python3

import os

import polars as pl
from genotations import ensembl
from genotations.quantification import *
import click
from genotations.genomes import Strand

"""
def search_in_expressions(df: pl.DataFrame,
                          gene_name: str,
                          tpm_columns: Union[pl.Expr, list[pl.Expr], "str", list[str]],
                          min_avg_value: float = 0.0,
                          exact: bool = True, genome: Optional[Genome] = None):
    search = pl.col("gene_name") == gene_name if exact else pl.col("gene_name").str.contains(gene_name)
    df = with_expressions_summaries(df, tpm_columns).filter(search).filter(pl.col("avg_TPM") >= min_avg_value)
    return df if genome is None else Annotations(df).with_sequences(genome).annotations_df
"""

def write_genes(bioprojects: Path, output: Path, skip_if_exists: bool = True, annotations: Optional[Annotations] = None,
                tpm_columns = pl.col("^SRR[a-zA-Z0-9]+$")):
    for p in bioprojects:
        name: str = p.stem + "_genes.parquet"
        where = output / name
        if where.exists() and skip_if_exists:
            print(f"{where} already exists, skipping...")
        else:
            expressions: pl.DataFrame = expressions_from_bioproject(p, False)
            print(f"results will be written to {where}")
            expressions.write_parquet(str(where))
            summarized_expressions =  with_expressions_summaries(expressions, tpm_columns)
            essential = [annotations.gene_col, annotations.gene_name_col]
            extended_expressions = Annotations(annotations.annotations_df.select(essential)).extend_with_annotations(summarized_expressions)
            extended_expressions.write_parquet(str(output / (p.stem + "_extended_genes.parquet")))


def write_transcript(p: Path, output: Path, skip_if_exists: bool = True,
                     annotations: Optional[Annotations] = None, genome: Optional[Genome] = None,
                     tpm_columns = pl.col("^SRR[a-zA-Z0-9]+$")):
    name = p.stem + "_transcripts.parquet"
    where = output / name
    expressions: pl.DataFrame = None
    if where.exists() and skip_if_exists:
        print(f"{where} already exists, skipping...")
        expressions = pl.read_parquet(where)
    else:
        expressions: pl.DataFrame = expressions_from_bioproject(p, True)
        print(f"results will be written to {where}")
        expressions.write_parquet(where)
    if annotations is not None:
        name_extended = p.stem + "_extended_transcripts.parquet"
        summarized_expressions = with_expressions_summaries(expressions, tpm_columns)
        essential = [annotations.gene_col, annotations.gene_name_col, annotations.transcript_col, annotations.transcript_name_col]
        extended_expressions = Annotations(annotations.annotations_df.select(essential)).extend_with_annotations(summarized_expressions).unique()
        extended_expressions.write_parquet(str(output / name_extended))
    if genome is not None:
        name_exons = p.stem + "_exons.parquet"
        with_exons_info = annotations.with_genes_transcripts_exons_coordinates_only().extend_with_annotations_and_sequences(summarized_expressions, genome, Strand.Undefined).drop("coordinates")
        with_exons_info.write_parquet(str(output / name_exons), use_pyarrow = True)
    print(output / name_exons, "was successfully written")


@click.command()
@click.option('--samples', type=click.Path(exists=True), help="samples folder", default="/data/samples/muscle_differentiation")
@click.option('--output', type=click.Path(), help="output")
@click.option('--skip_if_exist', type=click.BOOL, help="If it should skip if exists", default=True)
@click.option('--species', help="species", default = "mouse")
@click.option('--gtf', help="Optional path to GTF file")
def prepare(samples: str, output: str, skip_if_exist: bool = True,  species: str = "mouse", gtf: Optional[str] = None):
    bioprojects = dirs(Path(samples))
    from genotations.genomes import Annotations
    script_path = os.path.realpath(os.path.dirname(__file__))
    output_path = (Path(script_path) / ".." / "data" / "inputs").absolute().resolve() if output is None else Path(output).absolute().resolve()
    output_path.mkdir(exist_ok=True)
    annotations: Annotations = None
    genome: Optional[Genome] = None
    if gtf is not None:
        print(f"loading annotations from GTF {gtf}")
        annotations = Annotations(Path(gtf))
    if species is not None and species.lower() == "mouse" or species.lower() == "mus_musculus":
        print('importing annotations from mouse')
        annotations = ensembl.mouse.annotations
        genome = ensembl.mouse.genome
    elif species is not None and species.lower() == "human" or species.lower() == "homo_sapiens":
        print("importing annotations from human")
        annotations = ensembl.human.annotations
        genome = ensembl.human.genome
    else:
        raise Exception("No GTF and no species info provided, nothing to take annotations from!")
    print(f"saving results to {output_path}")
    write_genes(bioprojects, output_path, skip_if_exist, annotations)
    print("genes saved")
    for p in bioprojects:
        write_transcript(p, output_path, skip_if_exist, annotations, genome, tpm_columns = pl.col("^SRR[a-zA-Z0-9]+$"))
    print('transcripts saved')


if __name__ == '__main__':
    prepare()