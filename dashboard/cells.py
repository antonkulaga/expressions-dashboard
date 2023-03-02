#!/usr/bin/env python3
import dash
import pandas as pd
import polars as pl
from dash import dcc, html, Input, Output, dash_table
from pycomfort.files import *
from pprint import pprint
from dash.development.base_component import Component
from models import *
from genotations.genomes import Annotations
from genotations.primers import suggest_primers_for_transcript_by_exons, PrimerResults, PrimerResult
from config import Locations
from dash import Dash, dcc, html
from dash.dependencies import Input, Output, State, MATCH

app = dash.Dash(__name__, suppress_callback_exceptions=True)

base = Path("..") if (Path("..") / "dashboard").exists() else Path(".")
locations = Locations(base)

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
    dcc.Input(id='min_expression', type='number', value=0.1, min=0.1, max=1000000, step=0.1),
    html.H3('Gene expressions interface'),
    dcc.Tabs(id="tabs", value=bioprojects[0].id, children=seq(bioprojects).map(tab_from_bioproject).to_list()),
    html.Div(id='tabs-content')
])

def _add_links(row: list[str], k: any, v: str):
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

def df_to_table(df: pl.DataFrame, id: str = Component.UNDEFINED, of_type: str = "transcript_expressions") -> dash_table.DataTable:
    name_ids = [{"name": i, "id": i, "presentation": "markdown"} for i in df.columns]
    records: list[dict[str, any]] = df.to_dicts()
    for row in records:
        for k, v in row.items():
            _add_links(row, k, v)
            if isinstance(v, list):
                row[k] = str(v)
    return dash_table.DataTable(records, name_ids, id={
        "type": of_type,
        "index": id
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
                                },
                                filter_action="native",
                                sort_action="native",
                                sort_mode="multi",
                                column_selectable="single",
                                row_selectable="multi"
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
              Input('genes', 'value'),
              Input('min_expression', 'value')
              )
def render_content(tab: str, genes: list[str], min_avg_expression: float):
    assert tab in bioprojects_by_id, "bioproject and tab are connected"
    #print("GENES ARE: ", genes)
    bio: Bioproject = bioprojects_by_id[tab]
    runs_tab = df_to_table(bio.runs, tab+"_runs")
    transcripts_tab: dash_table.DataTable
    exons_tab: dash_table.DataTable
    selected_transcripts: pl. DataFrame
    available_exons: pl.DataFrame

    #exons_message = "Showing all available exons to choose from" if None else f"showing exons of selected transcript"

    if genes is None or len(genes) < 1:
        message = f"No genes selected, showing top 10 expressed genes instead"
        selected_transcripts = bio.transcripts.sort(pl.col("sum_TPM"), reverse=True).filter(pl.col("avg_TPM") >= min_avg_expression).head(10)
        available_exons = bio.exons.sort(pl.col("sum_TPM")).filter(pl.col("avg_TPM") >= min_avg_expression)
    else:
        message = f"No genes selected, showing top 10 expressed genes instead" if genes is None or len(genes) < 1 else f"Selected genes {genes}"
        selected_transcripts = bio.transcripts.filter(pl.col("gene_name").is_in(genes)).filter(pl.col("avg_TPM") >= min_avg_expression)
        available_exons = bio.exons.filter(pl.col("gene_name").is_in(genes)).filter(pl.col("avg_TPM") >= min_avg_expression)
    return html.Div([
        html.H3(f"{tab} content"),
        runs_tab,
        html.H4(message),
        df_to_heatmap(selected_transcripts, "Transcript expression heatmap", 40.0),
        html.H4("Transcript expression table"),
        df_to_table(selected_transcripts, tab+"_transcripts"),
        html.H4("Primer selection"),
        html.Plaintext("Transcript selection"),
        dcc.Dropdown(
            selected_transcripts.select(pl.col("transcript_name")).unique().to_series().to_list(),
            multi=False, id={
                "index": tab,
                "type": "transcript_for_primer"
            }, value=""
        ),
        df_to_table(available_exons, tab, of_type="selected_exons"),
        html.Plaintext("Exon for the forward primer"),
        dcc.Dropdown(
            available_exons.select(pl.col("exon_number")).unique().to_series().to_list(),
            multi=False, id={
                "index": tab,
                "type": "forward_exon"
            }, value=""
        ),
        html.Plaintext("Exon for the reverse primer"),
        dcc.Dropdown(
            available_exons.select(pl.col("exon_number")).unique().to_series().to_list(),
            multi=False, id={
                "index": tab,
                "type": "reverse_exon"
            }, value=""
        ),
        html.Plaintext("Optimal  temperature"), dcc.Input(id={
            "index": tab,
            "type": "optimal_temperature"
        }, type='number', value=60, min=58, max=65, step=0.5),
        html.Plaintext("Maximum amplified product"), dcc.Input(id={
            "index": tab,
            "type": "maximum_product"
        }, type='number', value=200, min=58, step=1),
        html.Br(),
        html.Button("Compute primers",
                    id={
                        "index": tab,
                        "type": "compute_button"
                    }, n_clicks=0
                    ),
        html.Br(),
        html.Output(id={
            "index": tab,
            "type": "computed_primers"
        })
    ])

def _strings_to_spans(strings: list[str]):
    return seq(strings).fold_left( ((0, 0),), lambda acc, el: acc + ((acc[-1][0]+acc[-1][1], len(el)),) ).to_list()[1:]

@app.callback(Output({'type': 'computed_primers', 'index': MATCH}, 'children'),
              State({'type': 'transcript_for_primer', 'index': MATCH}, 'value'),
              State({'type': 'forward_exon', 'index': MATCH}, 'value'),
              State({'type': 'reverse_exon', 'index': MATCH}, 'value'),
              State({'type': 'optimal_temperature', 'index': MATCH}, 'value'),
              State({'type': 'maximum_product', 'index': MATCH}, 'value'),
              State('tabs', 'value'),
              Input({'type': "compute_button", 'index': MATCH}, 'n_clicks')
              )
def compute_primers(transcript_for_primer: str,
                    forward_exon: int, reverse_exon: int,
                    optimal_temperature: float,
                    maximum_product: int,
                    tab: str, n_clicks: float):
    bio: Bioproject = bioprojects_by_id[tab]
    if transcript_for_primer is not None and reverse_exon > forward_exon:
        results: PrimerResults = suggest_primers_for_transcript_by_exons(Annotations(bio.exons), transcript_for_primer, forward_exon, reverse_exon, opt_t=optimal_temperature, max_product=maximum_product, debug=True)
        records = [{"forward_primer": r.PRIMER_LEFT_SEQUENCE,
                    "reverse_primer": r.PRIMER_RIGHT_SEQUENCE,
                    "left_primer_tm": r.PRIMER_LEFT_TM,
                    "right_primer_tm": r.PRIMER_RIGHT_TM,
                    "primer_product_size": r.PRIMER_PAIR_PRODUCT_SIZE
                    }  for r in results.get_cleaned_results()
                   ]
        return [dash_table.DataTable(records, id={
            "type": "primers_found",
            "index": tab
        })]
    else:
        return html.Plaintext("Please, select transcript for the primer!")


@app.callback(Output({'type': 'selected_exons', 'index': MATCH}, 'data'),
              Input({'type': 'transcript_for_primer', 'index': MATCH}, 'value'),
              Input('tabs', 'value'))
def update_exons(transcript_for_primer: str, tab: str):
    bio: Bioproject = bioprojects_by_id[tab]
    df = bio.exons.filter(pl.col("transcript_name") == transcript_for_primer)
    records: list[dict[str, any]] = df.to_dicts()

    for row in records:
        for k, v in row.items():
            _add_links(row, k, v)
    return records


if __name__ == '__main__':
    app.run_server(debug=True, host="0.0.0.0")