#!/usr/bin/env python3
import dash
from dash import Dash, dcc, html
from dash import dash_table
from dash.development.base_component import Component
from functional.pipeline import Sequence

from config import Locations
from cell_models import Cells
from models import *
from components import *
from dash import Dash, dcc, html
from dash.dependencies import Input, Output, State, MATCH, ALL
from functional import seq

app = dash.Dash(__name__, suppress_callback_exceptions=True)

base = Path("..") if (Path("..") / "dashboard").exists() else Path(".")
locations = Locations(base)
cells = Cells(locations.cell_lines)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    html.H1('Gene expressions interface'),
    df_to_table(cells.toc_samples, index="all", of_type="samples"),
    html.H1('Selected rows'),
    html.Div(id='selected_rows')
])

def _extract_samples(row: dict) -> Sequence:
    return seq(row["Samples"].split(",")).map(lambda v: v.strip())

@app.callback(
    Output('selected_rows', "children"),
    Input({'type': 'samples', 'index': "all"}, "derived_virtual_data"),
    Input({'type': 'samples', 'index': "all"}, "derived_virtual_selected_rows"))
def render(rows: dict, selected_row_indexes: dict):
    selected_indexes = [] if selected_row_indexes is None else selected_row_indexes
    if not selected_indexes:
        return html.Div(["Select cell lines!"])
    selected_rows = [rows[i] for i in selected_indexes]
    selected_samples = seq(selected_rows)\
        .flat_map(
            lambda r: _extract_samples(r).map(lambda s: {
                "sample": s,
                "cell_line": r["Cell line name"]
            }))\
        .to_list()
    selected_df = pl.from_dicts(selected_samples)
    table = df_to_table(selected_df, index="selected", of_type="samples")
    return table


if __name__ == '__main__':
    app.run_server(debug=True, host="0.0.0.0", port=8051)