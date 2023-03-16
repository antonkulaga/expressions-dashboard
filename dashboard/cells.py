#!/usr/bin/env python3
import genotations.ensembl
from dash import Dash, dcc, html
from dash.dependencies import Input, Output, MATCH, State
from functional import seq
from functional.pipeline import Sequence

from cell_models import Cells, Sample, SpeciesExpressions
from components import *
from config import Locations

app = dash.Dash(__name__, suppress_callback_exceptions=True)

base = Path("..") if (Path("..") / "dashboard").exists() else Path(".")
locations = Locations(base)
cells = Cells(locations.cell_lines, locations.cell_lines_cache)
external_scripts = [
    'https://cdn.jsdelivr.net/npm/jquery@3.6.3/dist/jquery.min.js',
    'https://cdn.jsdelivr.net/npm/fomantic-ui@2.9.2/dist/semantic.min.js'
]
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css',
                        "https://cdn.jsdelivr.net/npm/fomantic-ui@2.9.2/dist/semantic.min.css"
                        ]

heatmap_row = 40

app = Dash(__name__, external_stylesheets=external_stylesheets, external_scripts=external_scripts)

app.layout = html.Div([
    html.Div(
        [html.H1('Gene expressions interface')], className="ui top blue inverted segment"
    ),
    html.Div([
        "This user interface is devoted to exploring transcript expressions of cell lines that we have.",
        "For each cell line several samples were quantified. To view, select cell lines of your interest and then pick genes",
        "Expressions are output on per-species bases, you should open the accordion with corresponding species name"
    ], className="ui message"),
    html.Div([
    df_to_table(cells.toc_samples, index="cells", of_type="samples")]),
    html.Div([
        html.H1('Transcript expressions for selected runs'),
        html.Div([
            html.Div(["Minimal average expression in at least one run"], className="ui label"),
            dcc.Input(id='min_expression', type='number', value=0.1, min=0.1, max=1000000, step=0.1)
        ], className="ui labeled input"),
        html.Hr(className="ui divider"),
        html.Div(id='selection_results')
    ], className="ui blue segment")
])

def render_species_segment(species: SpeciesExpressions, selected_runs_df: pl.DataFrame):
    options = species.gene_names.select("gene_name").to_series().to_list()
    selected_runs_table = df_to_table(selected_runs_df, f"Selected runs with {species.species}", row_selectable=False)
    return [
        html.Summary([
            html.I(className="dropdown icon"),
           f'Sequencing runs for {species.species}'
        ], className="title"),
        html.Div(
            [ selected_runs_table,
              html.Div(["Select genes of your interest:"], className="ui label"),
              dcc.Dropdown(options, multi=True, id={ "type": "genes", "index": species.species}),
              html.H3("Transcript expressions:"),
              html.Div( id={ "type": "expressions", "index": species.species })
              ], className="content"
        ),
    ]



@app.callback(
    Output({'type': 'expressions', 'index': MATCH}, "children"),
    Input({'type': 'samples', 'index': "cells"}, "derived_virtual_data"),
    Input({'type': 'samples', 'index': "cells"}, "derived_virtual_selected_rows"),
    Input('min_expression', 'value'),
    Input({'type': "genes", "index": MATCH}, "value"),
    State({'type': "genes", "index": MATCH}, 'id'),
)
def render_species(rows: dict, selected_row_indexes: dict, min_avg_expression: float, genes: list[str], id: dict[str, str]):
    species_name = id["index"]
    if species_name is None:
        return html.Div(["Species is not selected!"])
    if species_name not in genotations.ensembl.species:
        return html.Div([f"Species {species_name} is not selected!"])
    no_genes = genes is None or len(genes) < 1
    print(f"GENES ARE {genes}")
    selected_species = cells.species_expressions_by_name[species_name]
    selected_indexes = [] if selected_row_indexes is None else selected_row_indexes
    selected_rows = [rows[i] for i in selected_indexes]
    selected_samples = seq(selected_rows).flat_map(lambda row: cells.extract_samples_from_row(row))
    selected_run_list: list[str] = selected_samples.flat_map(lambda s: [run.run_accession for run in s.run_list if run.scientific_name.replace(" ", "_") == species_name]).to_list()
    if len(selected_run_list) < 1:
        return html.Div(f"no runs for {species_name} found!")
    cols = seq(selected_species.transcript_expressions.columns)\
        .filter_not(lambda c: c.startswith("SRR") and c not in selected_run_list)\
        .to_list()
    transcripts_with_samples = selected_species.transcript_expressions.select(cols).filter(pl.col("avg_TPM") >= min_avg_expression)
    top_genes_number = 25
    selected_transcripts = transcripts_with_samples.head(top_genes_number).sort(pl.col("sum_TPM"), reverse=True) if no_genes \
        else transcripts_with_samples.filter(pl.col("gene_name").is_in(genes)).sort(pl.col("sum_TPM"), reverse=True)
    #preparing tables
    heatmap = df_to_heatmap(selected_transcripts, f"{species_name} samples expression heatmap", heatmap_row)
    species_table = df_to_table(selected_transcripts, selected_species.species+"_transcripts", row_selectable=False)
    message = f"as no genes were selected, showing top {top_genes_number} transcripts by avg expression" if no_genes else f"Expression of {', '.join(genes)}"
    return html.Div([
        html.H4(message),
        f"Only transcripts with minimal expressions higher than {min_avg_expression} TPM are shown",
        heatmap,
        species_table
                     ])

@app.callback(
    Output('selection_results', "children"),
    Input({'type': 'samples', 'index': "cells"}, "derived_virtual_data"),
    Input({'type': 'samples', 'index': "cells"}, "derived_virtual_selected_rows")
)
def render(rows: dict, selected_row_indexes: dict):
    selected_indexes = [] if selected_row_indexes is None else selected_row_indexes
    if not selected_indexes:
        return html.Div(["Select cell lines!"])
    selected_rows = [rows[i] for i in selected_indexes]
    selected_samples: list[Sample] = seq(selected_rows).flat_map(lambda row: cells.extract_samples_from_row(row)).to_list()
    species_names = set([r.scientific_name.replace(" ", "_") for s in selected_samples for r in s.run_list])
    species: list[SpeciesExpressions] = [cells.species_expressions_by_name[species_name] for species_name in species_names if species_name in cells.species_expressions_by_name]
    for s in species:
        print(f"SPECIES {s.species}")
    results = seq([
        render_species_segment(s, Sample.samples_to_df(selected_samples).filter(pl.col("scientific_name").str.contains(s.species))) for s in species if s.gene_expressions is not None and s.transcript_expressions is not None
    ]).flatten().to_list()
    header = html.H4(["Expressions by species:"], className="ui large header")
    return html.Div([
        header,
        html.Details(results, className="ui blue styled fluid accordion")
        ]
    )


if __name__ == '__main__':
    app.run_server(debug=True, host="0.0.0.0", port=8051)