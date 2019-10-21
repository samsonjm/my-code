#! /usr/bin/python3

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
import rpy2
import numpy as np
from sklearn.preprocessing import scale
import pandas as pd
from pandas import DataFrame
import re
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter

app = dash.Dash()
FILES = {'411 + 401':'syngenome.sleuth',
         '411 + 401 + Shoot':'basegenome.sleuth',
         '411':'411genome.sleuth',
         '411 + Shoot':'411shgenome.sleuth',
         '401':'401genome.sleuth',
         '401 + Shoot':'401shgenome.sleuth'}
readRDS = ro.r['readRDS']
MATRICES = {'syngenome.sleuth':readRDS('syngenome.sleuth'),
            'basegenome.sleuth':readRDS('basegenome.sleuth'),
            '411genome.sleuth':readRDS('411genome.sleuth'),
            '411shgenome.sleuth':readRDS('411shgenome.sleuth'),
            '401genome.sleuth':readRDS('401genome.sleuth'),
            '401shgenome.sleuth':readRDS('401shgenome.sleuth')}

app.layout = html.Div(children=[
    html.H1(children= '''Expression Series'''),
    dcc.Dropdown(options = [{'label': key, 'value': value} for key,
                             value in FILES.items()],
                 value='syngenome.sleuth',
                 id = "file"),
    dcc.Dropdown(id = "locus", multi=True, value=None),
    dcc.Checklist(id = "scale", options = [{'label':'Scale', 'value':'scaled'}]),
    dcc.Checklist(id = "corr", options = [{'label':'Correlation', 'value':'correlation'}]),
    dcc.Checklist(id = "absolute", options = [{'label':'Correlate Based on Absolute Value', 'value':'absol'}]),
    dcc.Slider(id = "cutoff", min=0, max=1, step=0.01, value=1,
                marks={0:'0', 0.25:'0.25', 0.50:'0.50', 0.75:'0.75', 1:'1'}),
    html.Div(children='''Output'''),
    dcc.Graph(id="expression")
    # html.Div(id='my-div')
])

@app.callback(
    Output('locus', 'options'),
    [Input('file', 'value')])
def populate_gene_names(input_file):
    if input_file == None:
        return None
    matrix = MATRICES[input_file]
    return [{'label': gene, 'value': gene} for gene in matrix.rownames]

@app.callback(
    Output('expression', 'figure'),
    [Input('file', 'value'),
     Input('locus', 'value'),
     Input('scale', 'value'),
     Input('corr', 'value'),
     Input('absolute', 'value'),
     Input('cutoff', 'value')])
def generate_expression_graph(input_file, loci, scale, correlation, absval, cutoff):
    if not loci:
        return {
            'data':[]
        }

    # Reads in the .sleuth RDS file
    matrix = MATRICES[input_file]

    # Converts the .RDS matrix into a dataframe
    with localconverter(ro.default_converter + pandas2ri.converter):
        pd_matrix = ro.conversion.rpy2py(matrix)
    df = pd.DataFrame(pd_matrix, index=matrix.rownames,
                      columns=matrix.colnames).T

    # Adds correlation functionality
    if(correlation == ['correlation'] and len(loci) == 1):
        expr_gene = df.T.loc[loci]
        if absval == ['absol']:
            genes = df.T.loc[np.abs(np.corrcoef(expr_gene, df.T)[0][1:]) >= cutoff].index
        else:
            genes = df.T.loc[np.corrcoef(expr_gene, df.T)[0][1:] >= cutoff].index
    else:
       genes = loci

    df = df[genes]
    construct = [re.sub("-.*", "", _) for _ in df.index]
    df['construct'] = construct
    time = [re.sub(".*_([0-9]+[hd]).*", "\\1", _) for _ in df.index]
    df['time'] = time
    treatment = [re.sub(".*_([a-zA-Z]*)[1-9]*", "\\1", _) for _ in df.index]
    df['treatment'] = treatment
    df['stage'] = df['treatment'] + " " + df['time']

    # Sorts the rows in the dataframe
    hourcounts = {"h":1, "d":24}
    treatmentcounts = {"Mock":0, "Mock ":0, "Dex":10000, "Dex ":10000, "Sh":1000000, "Sh ":1000000}
    sortvalue = [int(treatmentcounts[v[:-3]] + int(v[-3:-1])*hourcounts[v[-1]]) for v in df.stage]
    df['sortvalue'] = sortvalue
    df = df.sort_values('sortvalue')
    df = df.drop('sortvalue', axis=1)

    # Finishes preparing the dataframe for plotting
    df = df.melt(id_vars=["stage", "construct", "time", "treatment"])
    if (len(df.columns) == 1):
        return()
    df.columns = ["stage", "construct", "time", "treatment", "TAIR", "expression"]
    df['tairstage'] = df['TAIR'] + " " + df['construct'] + " " + df['stage']

    df['expmean'] = df.tairstage.map(df.groupby('tairstage').expression.mean())

    # Creates a non-smooth line
    df['zscore'] = df.groupby('TAIR')['expression'].apply(lambda x: (x - x.mean())/x.std())
    df['zmean'] = df.tairstage.map(df.groupby('tairstage').zscore.mean())

    # Adds scaling functionality
    mydata = []
    colors = ["red", "blue", "black", "orange", "yellow", "purple", "green", "brown"]

    if(scale == ['scaled']):
        scatteryvariable = 'zscore'
        lineyvariable = 'zmean'
    else:
        scatteryvariable = 'expression'
        lineyvariable = 'expmean'

    colorindex = 0
    for locus in genes:
        for inducer in df.construct.unique():
            color = colors[colorindex % 8]
            data = {'x':df[(df['TAIR']==locus) & (df['construct']==inducer)]['stage'].tolist(),
                'y':df[(df['TAIR']==locus) & (df['construct']==inducer)][scatteryvariable].tolist(),
                'type':'scattergl',
                'name':inducer + " " + locus,
                'mode':'markers',
                'marker':{'color':color}}
            mydata.append(data)
            data = {'x':df[(df['TAIR']==locus) & (df['construct']==inducer)]['stage'].tolist(),
                'y':df[(df['TAIR']==locus) & (df['construct']==inducer)][lineyvariable].tolist(),
                'type':'scattergl',
                'showlegend':False,
                'marker':{'color':color},
                'name':inducer + " " + locus}
            mydata.append(data)
            colorindex += 1

    figure={
        'data': mydata,
        'layout': {
            'title':input_file}
    }
    return figure

if __name__ == "__main__":
    app.run_server(debug=False, host='0.0.0.0', port = 54324)
