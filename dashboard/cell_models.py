from typing import *
from pathlib import Path

import genotations.ensembl
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
    def __init__(self, folder: Path, autoload_runs: bool = False):
        self.folder = folder.absolute().resolve()
        self.meta = pl.read_csv(folder / f"{folder.name}.tsv", sep="\t")
        self.parent_name = self.folder.parent.name

        self.run_list = [Run(d, folder, autoload_runs) for d in self.meta.to_dicts()]
        #self.runs = collections.OrderedDict([(run.run_accession, run) for run in self.run_list])


    def load(self) -> Run:
        for run in self.run_list:
            if run.genes is None:
                run.load()
        return self



class Cells:

    folder: Path
    toc_path: Path
    toc: pl.DataFrame
    toc_samples: pl.DataFrame

    def __init__(self, folder: Path):
        self.folder = folder.absolute().resolve()
        self.toc_path = self.folder / "cell_lines_toc.tsv"
        self.toc = pl.read_csv(str(self.toc_path), sep="\t")
        self.toc_samples = self.toc.filter(pl.col("Samples").str.contains("SAMN"))

    def extract_samples_from_row(self, row: dict, autoload: bool = False) -> Sequence:
        cell_line: str = row["Cell line name"]
        cell_line_folder = (self.folder / cell_line.strip()).absolute().resolve()
        if cell_line_folder.exists() == False:
            print(f"cell line {cell_line_folder} does not seem to exist!")
            return []
        return seq(row["Samples"].split(",")).map(lambda s: cell_line_folder / s)\
            .filter(lambda s: s.exists() and (s / f"{s.name}.tsv").exists())\
            .map(lambda f: Sample(f, autoload))

    def extract_all_samples(self, autoload: bool = False) -> list[Sample]:
        return seq(self.toc_samples.to_dicts()).flat_map(lambda s: self.extract_samples_from_row(s, autoload)).to_list()


class ExpressionSelection:

    selected_samples: list[Sample]
    runs_list: List[Run]
    runs_by_species: OrderedDict[str, List[Run]]
    gene_expressions: OrderedDict[str, pl.DataFrame]
    transcript_expressions: OrderedDict[str, pl.DataFrame]

    def __init__(self, selected_samples: list[Sample]):
        self.selected_samples = selected_samples
        self.runs_list = seq(self.selected_samples).flat_map(lambda s: s.run_list).to_list()
        self.runs_by_species = collections.OrderedDict(seq(self.runs_list).group_by(lambda r: r.scientific_name.replace(" ", "_")).to_list())
        acc_genes: list[(str, pl.DataFrame)] = [] #yes, it is ugly, I know
        acc_transcripts: list[(str, pl.DataFrame)] = [] #yes, it is ugly, I know
        for s, runs in self.runs_by_species.items():
            if s not in genotations.ensembl.species:
                print(f"{s} was not found in genotations!")
            else:
                species_info: genotations.ensembl.SpeciesInfo = genotations.ensembl.species[s]
                genes_for_merge = collections.OrderedDict(seq(runs).map(lambda r: (r.run_accession, r.load().genes)).to_list())
                merged_genes = genotations.quantification.merge_expressions(genes_for_merge, False)
                acc_genes.append((species_info.species_name, merged_genes))
                transcripts_for_merge = collections.OrderedDict(seq(runs).map(lambda r: (r.run_accession, r.load().transcripts)).to_list())
                merged_transcripts = genotations.quantification.merge_expressions(transcripts_for_merge, True)
                acc_transcripts.append((species_info.species_name, merged_transcripts))
        self.gene_expressions = collections.OrderedDict(acc_genes)
        self.transcript_expressions = collections.OrderedDict(acc_transcripts)




