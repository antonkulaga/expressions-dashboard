table_conditional_style = [
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
]

table_header_style = {
    'backgroundColor': 'rgb(210, 210, 210)',
    'color': 'black',
    'fontWeight': 'bold'
}