# import functions and text
from bin.dashboard.dashboard_functions import read_compsite_files, read_DFI_csvs, read_hng_files
from bin.dashboard.text_boxes import intro_md, flex_md, composite_md
import os

# File management/OS
from pathlib import Path, PurePosixPath


#Plotting
# import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

## DASH
from dash import Dash, dcc, html, Input, Output
import webbrowser
import subprocess


## Config

#external_stylesheets = ["./bWLwgP.css"]
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

all_outputs = []
for child in Path("test_output/").iterdir():
    if child.is_dir():
        all_outputs.append(str(child))

## Buttons, etc
app.layout = html.Div([
    html.Div(children=[
        html.Label('Select Output'),
        dcc.Dropdown(
            options=all_outputs,
            value=all_outputs[0],
            id="frag-dropdown",
            multi=False),

        dcc.Markdown(
            children=intro_md
        ),
        dcc.Graph(
            id='coverage-plot',
        ),
        
        html.Br(),
        html.Label('Radio Items'),
        ], 
    style={'padding': 10, 'flex': 1}
    ),
    
    html.Div(children=[
        dcc.Markdown(
            children=flex_md
        ),
        dcc.Graph(
            id='flex-plot',
        ),
    ], style={'padding': 10, 'flex': 1}),

    html.Div(children=[
        dcc.Markdown(
            children=composite_md
        ),
        dcc.Graph(
            id='composite-plot',
            # figure=fig3,
        ),
    ], style={'padding': 10, 'flex': 1}),

    html.Div(children=[
        # html.Label('Checkboxes'),
        # dcc.Checklist(
        #     options=[
        #         {"label" : x, "value" : x, "disabled" : False}
        #         for x in all_filenames
        #     ],
        #     value=[all_filenames[0]],
        #     id="frag-checkbox"
        # ),

        html.Br(),
        html.Label('Slider'),
        dcc.Slider(
            min=0,
            max=9,
            marks={i: f'Label {i}' if i == 1 else str(i) for i in range(1, 6)},
            value=5,
        )
    ], style={'padding': 10, 'flex': 1})])

@app.callback(
    Output('coverage-plot', 'figure'),
    [Input('frag-dropdown', 'value')]
    )
def update_graph(options_chosen):
    ## COVERAGE ##
    i = 0
    df_list = []
    structure_list = []
    for child in Path(os.path.join(options_chosen, "REPORT", "COVERAGE")).iterdir():
        if child.is_file() and "composite" not in str(child):
            i += 1
            df = pd.read_csv(child)
            df_list.append(df)
            structure_list.append(child)
                
    fig1 = make_subplots(rows=i, cols=1, shared_xaxes=True)

    i = 1
    for df in df_list:
        fig1.append_trace(go.Scatter(
            x=df[df.columns[0]], # ResID
            y=df[df.columns[1]],
            fill='tozeroy', 
            name=str(structure_list[i-1])
        ), row=i, col=1)
        i +=1

    fig1.update_layout(height=400, width=1000, title_text="Coverage")
    fig1.update_yaxes(showgrid=False, range=[0,1], nticks=2)
    
    return fig1





if __name__ == "__main__":
    # For Development only, otherwise use gunicorn or uwsgi to launch, e.g.
    # gunicorn -b 0.0.0.0:8050 index:app.server

    port = 8050 # Default port, change if occupied

    try:
        app.run_server(debug=True, port = port)
    
    except OSError:
        print(f""" Port already in use, in bin/make_dashboard.py line XX, 
    change the number of the port to one that is not occupied.

    If this is not the first time running this script, the port might not be 
    actually closed (maybe you used Ctrl+Z insteasd of Ctrl+C to exit the last 
    execution? ).

    Always terminate the execution of Dash with Ctrl+C
    ' """)
        try:
            print("Executing:  'fuser -k 8050/tcp'") 
            subprocess.run("fuser", "-k", f"{port}/tcp")
       
        except Exception:
            print("'fuser -k 8050/tcp' did not work, exiting")
            print(""" To check the ports already in use, try one of the following:
            
sudo lsof -i -P -n | grep LISTEN
sudo netstat -tulpn | grep LISTEN
sudo ss -tulpn | grep LISTEN
sudo lsof -i:22 ## see a specific port such as 22 ##
sudo nmap -sTU -O IP-address-Here

to free 127.0.0.1:8051 for example, type:

fuser -k 8051/tcp'

""")
            exit(1)

        
    
