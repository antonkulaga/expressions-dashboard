from typing import *
from pathlib import Path
import polars as pl
from pycomfort.files import *
from functional import seq

def deep_find(p: Path, fun: Callable[[Path], bool], depth: int = -1):
    fl = files(p).filter(fun)
    folds = dirs(p).filter(fun)
    inside = dirs(p).map(lambda d: deep_find(d, fun, depth - 1)) if depth != 0 else seq([])
    return fl.union(folds).union(inside)


class Cells:

    folder: Path
    toc_path: Path
    toc: pl.DataFrame
    toc_samples: pl.DataFrame

    def __init__(self, folder: Path):
        self.folder = folder
        self.toc_path = self.folder / "cell_lines_toc.tsv"
        self.toc = pl.read_csv(str(self.toc_path), sep="\t")
        self.toc_samples = self.toc.filter(pl.col("Samples").str.contains("SAMN"))

