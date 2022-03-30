from dashboard_functions import read_compsite_files, read_DFI_csvs, read_hng_files
from text_boxes import intro_md, flex_md, composite_md

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


## Cofig

#external_stylesheets = ["./bWLwgP.css"]
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
all_filenames = ["Exp1", "Exp2", "AF1", "AF2"]
app = Dash(__name__, external_stylesheets=external_stylesheets)
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True


### Data for the plots ###

## COVERAGE ##
i = 0
df_list = []
structure_list = []
for child in Path('/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/REPORT/COVERAGE').iterdir():
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


## HINGES/FLEXIBILITY ##
dfi_dict = read_DFI_csvs('/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/REPORT/DFI/')        
hng_dict = read_hng_files('/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/HINGES/')

fig2 = make_subplots(rows=len(dfi_dict.keys()), cols=1, shared_xaxes=True)

i = 1
for dfi_file in dfi_dict.keys():
    df = dfi_dict[dfi_file]
    fig2.append_trace(go.Scatter(
        x=df[df.columns[0]], # resIDs
        y=df[df.columns[1]], # pctdfi
        name=str(structure_list[i-1])
    ), row=i, col=1)
    if "AF_DFI" not in str(dfi_file.stem) and "RF_DFI" not in str(dfi_file.stem):
        j =1
        for hng_file in hng_dict.keys():
            if str(PurePosixPath(dfi_file).stem)[0:-13] == str(PurePosixPath(hng_file).stem):
                for hinge in hng_dict[hng_file]:
                    fig2.add_vrect(
                        x0=hinge.split(':')[0], 
                        x1=hinge.split(':')[1],
                        annotation_text=f"H{j}", annotation_position="top left",
                        fillcolor="#52BE80", opacity=0.2,
                        layer="below", line_width=0, 
                        row=i, col=1)
            j += 1
    i +=1
    fig2.update_layout(height=600, width=1200, title_text="DFI profiles + Predicted hinges", 
                      margin_pad=0, barmode="group")
    fig2.update_yaxes(showgrid=False, range=[0,1], nticks=2)

## COMPOSITE ##
composite_dir = "/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/REPORT/COVERAGE/"
comp_dict = read_compsite_files(composite_dir)


# Plot
comp_filenames = []
fig3 = make_subplots(rows=len(comp_dict.keys())+1, cols=1, shared_xaxes=True)

for file in comp_dict.keys():
    df = comp_dict[file]
    print(len(df.columns))
    fig3 = make_subplots(rows=len(df.columns)-1, cols=1, shared_xaxes=True)
    i = 0
    for column in df.columns:
        if i >= 1:
            fig3.append_trace(go.Scatter(
                x=df.iloc[:,0],
                y=df[df.columns[i]],
                fill='tozeroy',
                name=str(column)
            ), row=i, col=1)
        i +=1
        comp_filenames.append(file)
    fig3.update_layout(height=600, width=1200, title_text="Composite coverage", 
                          margin_pad=0, barmode="overlay")
    fig3.update_yaxes(showgrid=False, range=[0,1], nticks=2)

    
all_filenames = []
for df in df_list:
    i = 0
    for column in df.columns:
        if i >= 1:
            all_filenames.append(str(column)) 
        i+=1





## Buttons, etc
app.layout = html.Div([
    html.Div(children=[
        html.Label('Select Output'),
        dcc.Dropdown(
            all_filenames,
            id="frag-dropdown",
            multi=True),

        dcc.Markdown(
            children=intro_md
        ),
        dcc.Graph(
            id='coverage-plot',
            figure=fig1,
        ),
        
        
        html.Br(),
        html.Label('Radio Items'),
        dcc.RadioItems(
            ['New York City', 'Montréal', 'San Francisco'], 'Montréal'
            ),
    ], style={'padding': 10, 'flex': 1}),
    
    html.Div(children=[
        dcc.Markdown(
            children=flex_md
        ),
        dcc.Graph(
            id='flex-plot',
            figure=fig2
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
        html.Label('Checkboxes'),
        dcc.Checklist(
            options=[
                {"label" : x, "value" : x, "disabled" : False}
                for x in all_filenames
            ],
            value=[all_filenames[0]],
            id="frag-checkbox"
        ),

        html.Br(),
        html.Label('Slider'),
        dcc.Slider(
            min=0,
            max=9,
            marks={i: f'Label {i}' if i == 1 else str(i) for i in range(1, 6)},
            value=5,
        )
    ], style={'padding': 10, 'flex': 1})])


# @app.callback(
#     Output('composite-plot', 'figure'),
#     [Input('frag-checkbox', 'value')]
#     )
# def update_figure(values):
#     fig3 = make_subplots(rows=len(comp_dict.keys())+1, cols=1, shared_xaxes=True)
#     for file in comp_dict.keys():
#         df = comp_dict[file]
#         fig3 = make_subplots(rows=len(df.columns)-1, cols=1, shared_xaxes=True)
#         i = 0
#         for column in df.columns:
#             if i >= 1:
#                 fig3.append_trace(go.Scatter(
#                     x=df.iloc[:,0],
#                     y=df[df.columns[i]],
#                     fill='tozeroy',
#                     name=str(column)
#                 ), row=i, col=1)
#             i +=1
#         fig3.update_layout(height=600, width=1200, title_text="Composite coverage", 
#                             margin_pad=0, barmode="overlay")
#         fig3.update_yaxes(showgrid=False, range=[0,1], nticks=2)
    

#     fig3.update_layout(transition_duration=500)

#     return fig3


@app.callback(
    Output('composite-plot', 'figure'),
    [Input('frag-checkbox', 'value')]
    )
def update_graph(options_chosen):
    # for file in options_chosen:
    fig3 = make_subplots(rows=len(options_chosen), cols=1, shared_xaxes=True)
    for file in options_chosen:
        df = comp_dict[file]
        i = 1
        for column in df.columns:
            fig3.append_trace(go.Scatter(
                x=df.iloc[:,0],
                y=df[df.columns[i]],
                fill='tozeroy',
                name=str(column)
            ), row=i, col=1)
            i +=1
        fig3.update_layout(height=600, width=1200, title_text="Composite coverage", 
                            margin_pad=0, barmode="overlay")
        fig3.update_yaxes(showgrid=False, range=[0,1], nticks=2)
    return fig3



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

        
    
