import dataclasses
from typing import *
from pathlib import Path

import genotations.ensembl
import pandas as pd
import polars as pl
from pycomfort.files import *
from functional import seq
from genotations.quantification import read_quant
import collections


def quant_from_run(run: Path, name_part: str = "quant.sf", dir_part: str = "quant_") -> Optional[pl.DataFrame]:
    transcripts = "gene" not in dir_part and "gene" not in name_part
    qq = dirs(run).filter(lambda f: dir_part in f.name)
    if qq.len() < 1:
        print(f"could not find quantification data for {run}")
        return None
    else:
        if qq.len() > 1:
            print("there are more than two quant files, we are taking the first one")
        q = qq.first()
        return read_quant(files(q).filter(lambda f: name_part in f.name).first(), transcripts)


class Run:
    study_accession: str
    secondary_study_accession: str
    sample_accession: str
    secondary_sample_accession: str
    experiment_accession: str
    run_accession: str
    tax_id: str
    scientific_name: str
    instrument_model: str
    library_name: str
    library_layout: str
    library_strategy: str
    library_selection: str
    read_count: str
    experiment_title: str
    study_title: str
    experiment_alias: str
    fastq_ftp: str
    submitted_ftp: str
    sra_ftp: str
    sample_title: str
    folder: Path
    transcripts: pl.DataFrame
    genes: pl.DataFrame

    def __init__(self, dictionary: Dict[str, str], sample_folder: Path, autoload: bool = True):
        for k, v in dictionary.items():
            setattr(self, k, v)
        self.folder = (sample_folder / self.run_accession).absolute().resolve()
        assert self.folder, f"run folder {self.folder} does not exist!"
        self.quant_folder = self.folder / f"quant_{self.run_accession}"
        if autoload:
            self.load()

    def load(self) -> 'Run':
        self.transcripts = quant_from_run(self.folder)
        self.genes = quant_from_run(self.folder, name_part="quant.genes.sf")
        return self

class Sample:
    folder: Path
    meta: pl.DataFrame
    run_list: list[Run]
    runs_by_id: OrderedDict[str, Run]
    gene_expressions: pl.DataFrame
    transcript_expressions: pl.DataFrame
    parent_name: str

    @staticmethod
    def samples_to_df(samples: list['Sample']) -> Optional[pl.DataFrame]:
        metas = seq(samples).map(lambda s: s.meta).to_list()
        if len(metas) == 0:
            return None
        elif len(metas) == 1:
            return metas[0]
        else:
            return pl.concat(metas).with_column(pl.col("scientific_name").str.replace(" ", "_"))

    def __init__(self, folder: Path, autoload_runs: bool = False):
        self.folder = folder.absolute().resolve()
        self.meta = pl.read_csv(folder / f"{folder.name}.tsv", sep="\t")
        self.parent_name = self.folder.parent.name
        self.run_list = [Run(d, folder, autoload_runs) for d in self.meta.to_dicts()]


    def load(self) -> Run:
        for run in self.run_list:
            if run.genes is None:
                run.load()
        return self


class SpeciesExpressions:

    species: str
    gene_expressions: Optional[pl.DataFrame]
    transcript_expressions: Optional[pl.DataFrame]
    transcript_expressions_with_exons: Optional[pl.DataFrame]
    gene_names: Optional[pl.DataFrame]

    @staticmethod
    def load_summary_from(folder: Path, species: str):
        gene_expressions_path = folder / f"{species}_gene_expressions.parquet"
        transcript_expressions_path = folder / f"{species}_transcript_expressions.parquet"
        transcript_expressions_with_exons_path = folder / f"{species}_transcript_expressions_with_exons.parquet"
        gene_names_path = folder / f"{species}_gene_names.parquet"
        assert gene_expressions_path.exists(), f"gene expressions should exist for {species}"
        assert transcript_expressions_path.exists(), f"transcript expressions should exist for {species}"
        assert transcript_expressions_with_exons_path.exists(), f"{str(transcript_expressions_with_exons_path)} should exist for {species}"
        gene_names_df = pl.read_parquet(gene_names_path) if gene_expressions_path.exists() else None
        return SpeciesExpressions(species,
                           pl.read_parquet(str(gene_expressions_path)),
                           pl.read_parquet(str(transcript_expressions_path)),
                           pl.read_parquet(str(transcript_expressions_with_exons_path)),
                           gene_names_df
        )

    @staticmethod
    def load_summaries_from(folder: Path):
        species = files(folder).filter(lambda d: "_transcript_expressions.parquet" in d.name)\
            .map(lambda v: v.name.replace("_transcript_expressions.parquet", ""))
        return species.map(lambda s: SpeciesExpressions.load_summary_from(folder, s))

    def __init__(self, species: str,
                 gene_expressions: Optional[pl.DataFrame],
                 transcript_expressions: Optional[pl.DataFrame],
                 transcript_expressions_with_exons: Optional[pl.DataFrame],
                 gene_names: Optional[pl.DataFrame] = None
                 ):
        self.species = species.replace(" ", "_")
        self.gene_expressions = gene_expressions
        self.transcript_expressions = transcript_expressions
        self.transcript_expressions_with_exons = transcript_expressions_with_exons
        self.gene_names = self.annotator().gene_names_df.filter(pl.col("gene_name").is_not_null()) if gene_names is None else gene_names

    def write(self, folder: Path):
        folder.mkdir(exist_ok=True, parents=True)
        gene_expressions_path = folder / f"{self.species}_gene_expressions.parquet"
        transcript_expressions_path = folder / f"{self.species}_transcript_expressions.parquet"
        transcript_expressions_with_exons_path = folder / f"{self.species}_transcript_expressions_with_exons.parquet"
        gene_names_path = folder / f"{self.species}_gene_names.parquet"
        self.gene_expressions.write_parquet(str(gene_expressions_path), compression_level=22, use_pyarrow=True)
        self.transcript_expressions.write_parquet(str(transcript_expressions_path), compression_level=22, use_pyarrow=True)
        self.transcript_expressions_with_exons.write_parquet(str(transcript_expressions_with_exons_path), compression_level=22, use_pyarrow=True)
        self.annotator().gene_names_df.filter(pl.col("gene_name").is_not_null()).write_parquet(gene_names_path, compression_level=22, use_pyarrow=True)
        print(f"writing {gene_expressions_path} and {transcript_expressions_path}")
        return folder

    def annotator(self) -> genotations.genomes.Annotations:
        assert self.species in genotations.ensembl.species, f"species {self.species} should be in the genotations ensembl list!"
        return genotations.ensembl.species[self.species].annotations



    def is_empty(self):
        return self.transcript_expressions is None


    @staticmethod
    def empty(species: str):
        return SpeciesExpressions(species, None, None, None)


    @staticmethod
    def summaries_from_samples(selected_samples: list[Sample]):
        runs_list: list[Run] = seq(selected_samples).flat_map(lambda s: s.run_list).to_list()
        return SpeciesExpressions.summaries_from_runs(runs_list)

    @staticmethod
    def summaries_from_runs(runs_list: List[Run]):
        return seq(runs_list).group_by(lambda run: run.scientific_name.replace(" ", "_")) \
            .map(lambda v: SpeciesExpressions.from_runs(v[0], v[1])).to_list()


    @staticmethod
    def from_runs(species: str, runs: List[Run]) -> 'SpeciesExpressions':
        """
        Summarizes expressions from runs on species basis
        :param species:
        :param runs:
        :return:
        """
        assert species in genotations.ensembl.species, f"species f{species} was not found in genotations ensembl species dictionary!"
        species_info: genotations.ensembl.SpeciesInfo = genotations.ensembl.species[species]
        for_dic = seq(runs).map(lambda r: (r.run_accession, r.load().genes)).filter(lambda rg: rg[1] is not None)
        if for_dic.len() > 0:
            genes_for_merge = collections.OrderedDict(for_dic)
            merged_genes = genotations.quantification.merge_expressions(genes_for_merge, False)
            transcripts_for_merge = collections.OrderedDict(seq(runs).map(lambda r: (r.run_accession, r.load().transcripts)).to_list())
            merged_transcripts = genotations.quantification.merge_expressions(transcripts_for_merge, True)
            merged_summarized_transcripts =  genotations.quantification.with_expressions_summaries(merged_transcripts)
            merged_transcripts_extended = species_info.annotations.with_genes_transcripts_coordinates_only().extend_with_annotations(merged_summarized_transcripts)
            merged_transcripts_extended_exons = species_info.annotations.with_genes_transcripts_exons_coordinates_only().extend_with_annotations(merged_summarized_transcripts)
            return SpeciesExpressions(species, merged_genes, merged_transcripts_extended, merged_transcripts_extended_exons )
        else:
            try:
                print(f"no quantified transcripts found in {seq(runs).reduce(lambda a,b: a.run_accession+';'+b.run_accession)}")
            except:
                print(f"no quantified transcripts (second try) found in {seq(runs)}")
            return SpeciesExpressions.empty(species)


class Cells:

    folder: Path
    toc_path: Path
    toc: pl.DataFrame
    toc_samples: pl.DataFrame
    cache_folder: Path

    species_expressions: list[SpeciesExpressions]
    species_expressions_by_name: dict[str, SpeciesExpressions]


    def __init__(self, folder: Path):
        self.folder = folder.absolute().resolve()
        self.toc_path = self.folder / "cell_lines_toc.tsv"
        self.toc = pl.read_csv(str(self.toc_path), sep="\t")
        self.toc_samples = self.toc.filter(pl.col("Samples").str.contains("SAMN"))
        self.samples = self.extract_all_samples()
        self.cache_folder = folder / "CACHE"
        if self.cache_folder.exists() and files(self.cache_folder).len() >= 3:
            self.species_expressions = SpeciesExpressions.load_summaries_from(self.cache_folder)
        else:
            print(f"no cache folder found, producing expressions from samples and writing cache to {self.cache_folder}!")
            self.species_expressions = SpeciesExpressions.summaries_from_samples(self.samples)
            self.cache_folder.mkdir(exist_ok=True, parents=True)
            for s in self.species_expressions:
                if s.gene_expressions is None:
                    print(f"ERROR, {s.species} has None in gene expressions!")
                else:
                    s.write(self.cache_folder)
        self.species_expressions_by_name = {s.species: s for s in self.species_expressions}


    def extract_samples_from_row(self, row: dict, autoload: bool = False) -> Sequence:
        cell_line: str = row["Cell line name"]
        cell_line_folder = (self.folder / cell_line.strip()).absolute().resolve()
        if not cell_line_folder.exists():
            print(f"cell line {cell_line_folder} does not seem to exist!")
            return []
        return seq(row["Samples"].split(",")).map(lambda s: cell_line_folder / s.strip())\
            .filter(lambda s: s.exists() and (s / f"{s.name}.tsv").exists())\
            .map(lambda f: Sample(f, autoload))

    def extract_all_samples(self, autoload: bool = False) -> list[Sample]:
        return seq(self.toc_samples.to_dicts()).flat_map(lambda s: self.extract_samples_from_row(s, autoload)).to_list()







