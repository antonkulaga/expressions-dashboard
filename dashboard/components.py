from typing import Union

import polars as pl
from pathlib import *
from styles import *
import dash
from dash import Dash, dcc, html
from dash import dash_table
from dash.development.base_component import Component
from utils import *

def df_to_table(df: pl.DataFrame, index: str,
                of_type: str = "transcript_expressions", row_selectable: Union[str, bool] = "multi") -> dash_table.DataTable:
    name_ids = [{"name": i,
                 "id": i,
                 "presentation": "markdown"} for i in df.columns]
    records: list[dict[str, any]] = df.to_dicts()
    for row in records:
        for k, v in row.items():
            row[k] = make_clickable(k, v)
            if isinstance(v, list):
                row[k] = str(v)
    return dash_table.DataTable(records, name_ids, id={
        "type": of_type,
        "index": index
    },
        style_data_conditional=table_conditional_style,
        style_header=table_header_style,
        filter_action="native",
        sort_action="native",
        sort_mode="multi",
        column_selectable="single",
        row_selectable=row_selectable,
        style_table={'overflowX': 'scroll'}
        )

def df_to_heatmap(df: pl.DataFrame,  title: str = "Heatmap",
                  row_height: float = 40,
                  col: pl.Expr = pl.col("transcript_name"),
                  tpms: pl.Expr = pl.col("^SRR[a-zA-Z0-9]+$"),
                  ) -> dcc.Graph:
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
                    height= row_height  * (df.shape[0] +1), #cannot get it right for some reason
                    title=title
                    )
    fig.update_xaxes(side="top", tickfont = dict(size=20))
    return dcc.Graph(figure=fig, id = {
        "id": "heatmap",
        "index": title
    })