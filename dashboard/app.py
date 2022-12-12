#!/usr/bin/env python3
import dash
import pandas as pd
import polars as pl
from dash import dcc, html, Input, Output, dash_table
from pycomfort.files import *
from pprint import pprint
from dash.development.base_component import Component
from models import *


app = dash.Dash(__name__, suppress_callback_exceptions=True)

base = Path("..") if (Path("..") / "dashboard").exists() else Path(".")
data = base / "data"
inputs = data / "inputs"

bioprojects: list[Bioproject] = with_ext(inputs, ".tsv").map(lambda p: Bioproject(p.stem, inputs)).to_list()
bioprojects_by_id = {b.id: b for b in bioprojects}


from dash import Dash, dcc, html
from dash.dependencies import Input, Output

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = Dash(__name__, external_stylesheets=external_stylesheets)

def tab_from_bioproject(b: Bioproject):
    return dcc.Tab(label=b.id, value=b.id)

app.layout = html.Div([
    html.H1('Gene expressions interface'),
    dcc.Dropdown(
        bioprojects[0].gene_symbol_list,
        multi=True, id="genes", value=['Zswim8', 'Nf2', 'Dtx2', 'Traf3']
    ),
    dcc.Tabs(id="tabs", value=bioprojects[0].id, children=seq(bioprojects).map(tab_from_bioproject).to_list()),
    html.Div(id='tabs-content')
])

def df_to_table(df: pl.DataFrame, id: str = Component.UNDEFINED) -> dash_table.DataTable:
    name_ids = [{"name": i, "id": i, "presentation": "markdown"} for i in df.columns]
    records: list[dict[str, any]] = df.to_dicts()
    for row in records:
        for k, v in row.items():
            if k == "sample_accession":
                row[k] = f"[{v}](https://www.ncbi.nlm.nih.gov/biosample/{v})"
            elif k == "study_accession":
                row[k] = f"[{v}](https://www.ncbi.nlm.nih.gov/bioproject/?term={v})"
            elif k == "run_accession":
                row[k] = f"[{v}](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc={v})"
            elif k == "experiment_accession":
                row[k] = f"[{v}](https://www.ncbi.nlm.nih.gov/sra/SRX5864477/{v})"
            elif k == "gene":
                row[k] = f"[{v}](https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g={v})"
            elif k == "transcript":
                row[k] = f"[{v}](https://www.ensembl.org/Mus_musculus/Transcript/Summary?t={v})"
            elif k == "gene_name":
                row[k] = f"[{v}](https://www.genecards.org/cgi-bin/carddisp.pl?gene={v})"
    return dash_table.DataTable(records, name_ids, id={
        "type": "transcript_expressions",
        "id": id
    },
    style_data_conditional=[
        {
            'if': {'row_index': 'odd'},
            'backgroundColor': 'rgb(220, 220, 220)',
        },
        {
           'if': {'column_id': "sum_TPM"},
           'backgroundColor': 'rgb(255, 255, 204)',
        },
        {
            'if': {'column_id': "avg_TPM"},
            'backgroundColor': 'rgb(255, 255, 204)',
        }
    ],
    style_header={
        'backgroundColor': 'rgb(210, 210, 210)',
        'color': 'black',
        'fontWeight': 'bold'
    }
)

def df_to_heatmap(df: pl.DataFrame,  title: str = "Heatmap", row_height: float = 40, col: pl.Expr = pl.col("transcript_name"), tpms: pl.Expr = pl.col("^SRR[a-zA-Z0-9]+$")) -> dcc.Graph:
    from dash import dcc
    import plotly.express as px

    df_pandas = df.select(tpms).to_pandas()
    y = df.select(col).to_series().to_list()
    z_max = float(df_pandas.to_numpy().max())

    fig = px.imshow(df_pandas,
                    y = y,
                    zmin=0.0,
                    zmax=z_max,
                    text_auto=True,
                    aspect="auto",
                    height=row_height * df.shape[0],
                    title=title
                    )

    return dcc.Graph(figure=fig)


@app.callback(Output('tabs-content', 'children'),
              Input('tabs', 'value'),
              Input('genes', 'value')
              )


def render_content(tab: str, genes: list[str]):
    assert tab in bioprojects_by_id, "bioproject and tab are connected"
    #print("GENES ARE: ", genes)
    bio: Bioproject = bioprojects_by_id[tab]
    runs_tab = df_to_table(bio.runs, tab+"_runs")
    transcripts_tab: dash_table.DataTable
    exons_tab: dash_table.DataTable
    selected_transcripts: pl.DataFrame
    selected_exons: pl.DataFrame

    if genes is None or len(genes) < 1:
        message = f"No genes selected, showing top 10 expressed genes instead"
        selected_transcripts = bio.transcripts.sort(pl.col("sum_TPM"), reverse=True).head(10)
        selected_exons = bio.exons.sort(pl.col("sum_TPM"))
    else:
        message = f"No genes selected, showing top 10 expressed genes instead" if genes is None or len(genes) < 1 else f"Selected genes {genes}"
        selected_transcripts = bio.transcripts.filter(pl.col("gene_name").is_in(genes))
        selected_exons = bio.exons.filter(pl.col("gene_name").is_in(genes))
    return html.Div([
        html.H3(f"{tab} content"),
        runs_tab,
        html.H4(message),
        df_to_heatmap(selected_transcripts, "Transcript expression heatmap", 40.0),
        df_to_table(selected_transcripts, tab+"_transcripts")
    ])


if __name__ == '__main__':
    app.run_server(debug=True)