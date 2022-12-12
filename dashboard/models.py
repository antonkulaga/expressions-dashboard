import polars as pl
from functools import cache, cached_property
from pathlib import Path


class Bioproject:
    """
    Bioproject info for the UI
    """

    id: str
    folder: Path

    def __init__(self, id: str, folder: Path):
        self.id = id
        self.folder = folder

    @cached_property
    def runs(self) -> pl.DataFrame:
        p = self.folder / (self.id + ".tsv")
        assert p.exists(), f"tsv path {p} does not exist!"
        return pl.read_csv(p, sep="\t").select([
            "sample_accession",
            "sample_title",
            "run_accession",
            "study_accession",
            "experiment_alias",
            "experiment_accession",
            "scientific_name",
            "instrument_model",
            "library_name",
            "library_strategy",
            "library_selection",
            "read_count",
            "fastq_ftp"
        ])


    @cached_property
    def transcripts(self) -> pl.DataFrame:
        p = self.folder / (self.id + "_extended_transcripts.parquet")
        assert p.exists(), f"transcripts path {p} does not exist!"
        return pl.read_parquet(p)

    @cached_property
    def genes(self) -> pl.DataFrame:
        return self.transcripts.select([pl.col("gene"), pl.col("gene_name")]).unique()
    @cached_property
    def gene_symbol_list(self):
        return self.genes.select(pl.col("gene_name")).to_series().to_list()
    @cached_property
    def exons(self) -> pl.DataFrame:
        p = self.folder / (self.id + "_exons.parquet")
        assert p.exists(), "exons {p} does not exist!"
        return pl.read_parquet(p)