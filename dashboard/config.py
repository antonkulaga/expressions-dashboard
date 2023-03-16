from typing import *
from pathlib import Path
import polars as pl

def apply_settings():
    pl.Config.set_tbl_width_chars(10000)
    pl.Config.set_fmt_str_lengths(1000)
    pl.Config.set_tbl_rows(20)

class Locations:
    base: Path
    inputs: Path
    interim: Path
    output: Path
    cell_lines: Path
    cell_lines_cache: Path
    muscle_differentiation: Path
    def __init__(self, base: Path):
        self.base = base.absolute().resolve()
        self.data = base / "data"
        self.inputs = self.data / "inputs"
        self.inputs.mkdir(exist_ok=True)
        self.interim = self.data / "interim"
        self.interim.mkdir(exist_ok=True)
        self.output = self.data / "output"
        self.output.mkdir(exist_ok=True)
        self.cell_lines = self.inputs / "CELL_LINES"
        self.cell_lines.mkdir(exist_ok=True)
        self.muscle_differentiation = self.inputs / "muscle_differentiation"
        self.cell_lines_cache = self.inputs / "CACHE"
        self.cell_lines_cache.mkdir(exist_ok=True, parents=True)