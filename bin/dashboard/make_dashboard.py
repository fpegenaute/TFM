#Import all the needed python modules:

# Data Wrangling
import numpy as np
import pandas as pd

# File management/OS
import fnmatch
import os
from pathlib import Path, PurePosixPath

#Plotting
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go

## DASH
import dash
from dash import dcc
from dash import html
import webbrowser
import subprocess

# A couple of function definitions to make things easier later
def read_hng_files(hng_dir):
    """
    Given a folder containing PACKMAN .hng files, return a dictionary with the format:
    {filename : start:end residues}
    """
    i = 0
    hinges_dict = {}
    for child in Path(hng_dir).iterdir():
        if child.is_file() and fnmatch.fnmatch(child, "*.hng"):
            i += 1
            filename = os.path.basename(child)
            hinges_domains_df = pd.read_csv(child, sep="\t", names=["Chain", "Classification", "Start:End"])        
            hinge_df = hinges_domains_df[hinges_domains_df['Classification'].str.match('^H.*')== True]
            hinges_dict.update({filename[0:6] : hinge_df["Start:End"].tolist()})
    return hinges_dict

def read_DFI_csvs(dfi_csv_dir):
    """
    Given a folder containing csv files containing Dynamic Flexibility Index info, 
    return a dictionary with the format:
    {filename : 'start:end' residues}
    
    Format of the csv:
        header: Chain,ResID,pctdfi
    Chain example:3a58_A 
    """
    i = 0
    df_dict = {}
    for child in Path(dfi_csv_dir).iterdir():
        if child.is_file():
            i += 1
            df = pd.read_csv(child)
            df_dict.update({os.path.basename(child)[0:6] : df})
    return df_dict

def read_compsite_files(composite_dir):
    """
    Given a folder containing composite .csv files, return a dictionary with the format:
    {psition : coverage (0/1)}
    """
    i = 0
    comp_dict = {}
    for child in Path(composite_dir).iterdir():
        if child.is_file() and fnmatch.fnmatch(child, "*composite_coverage.csv"):
            i += 1
            filename = os.path.basename(child)
            df = pd.read_csv(child)
            comp_dict.update({os.path.basename(child)[0:6] : df})
    return comp_dict

## Cofig

#external_stylesheets = ["./bWLwgP.css"]
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
all_filenames = ["Exp1", "Exp2", "AF1", "AF2"]
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
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
        structure_list.append(os.path.basename(child))
            
fig1 = make_subplots(rows=i, cols=1, shared_xaxes=True)

i = 1
for df in df_list:
    fig1.append_trace(go.Scatter(
        x=df[df.columns[0]], # ResID
        y=df[df.columns[1]],
        fill='tozeroy', 
        name=structure_list[i-1]
    ), row=i, col=1)
    i +=1

fig1.update_layout(height=400, width=1000, title_text="Coverage")
fig1.update_yaxes(showgrid=False, range=[0,1], nticks=2)


## HINGES/FLEXIBILITY ##
dfi_dict = read_DFI_csvs('/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/REPORT/DFI/')        
hng_dict = read_hng_files('/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/HINGES/')

fig2 = make_subplots(rows=len(dfi_dict.keys()), cols=1, shared_xaxes=True)

i = 1
for file in dfi_dict.keys():
    df = dfi_dict[file]
    fig2.append_trace(go.Scatter(
        x=df["ResID"],
        y=df["pctdfi"],
        name=structure_list[i-1]
    ), row=i, col=1)
    j =1
    for hinge in hng_dict[file]:
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

# Read the files
composite_dir = "/home/gallegolab/Desktop/TFM/TFM/test_output/SEC3/REPORT/COVERAGE/"
comp_dict = read_compsite_files(composite_dir)

# Plot
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


# Text boxes
intro_md = '''
# Report

This is a summary of the results 

## Coverage

In this grapph,you can see whic parts of the reference FASTA sequence are 
covered by structures in the PDB
'''

flex_md = """
## Hinges and flexibility

Flexibility is an important feature of proteins, since they need to move to 
perform their function and interact with their substrates. In the following 
section, we provide you with two types of flexibility prediction: the Dynamic 
Flexibility Index and Hinge Prediction.

*Dynamic Flexibility Index*  
This is per-residue index indicating the contribution of each residue to the 
overall flexibility of the protein. It uses a method based in an Elastic Network 
Model, which is a more lightweight (but less precise, obviously) alternative to 
Molecular Dynamics. for ore info, 
[here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3673471/) is the original 
paper.

*Hinge Prediction*  
Hinges are the regions of the protein that allow it to move and change 
conformations. Using 
[this tool](https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbac007/6525212?login=true) 
the predicted hinge regions are showed on top of the DFI plot, with the 
significative ones colored in green, and  the non-significative ones in red.

"""

composite_md = """
## Composite

From all the structures retrieved by the program and provided by the user, 
make a composite with all of them that covers as much of the reference sequence 
as possible (avoiding overlaps).

This comopsite will be used to automatically build an 
[IMP topology file](https://integrativemodeling.org/2.5.0/doc/ref/classIMP_1_1pmi_1_1topology_1_1TopologyReader.html)
"""


## Buttons, etc
app.layout = html.Div([
    html.Div(children=[
        dcc.Markdown(
            children=intro_md
        ),
        dcc.Graph(
            id='coverage-plot',
            figure=fig1
        ),
        
        html.Label('Multi-Select Custom fragments'),
        dcc.Dropdown(
            all_filenames,
            ['Montréal', 'San Francisco'],
            multi=True),

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
            figure=fig3
        ),
    ], style={'padding': 10, 'flex': 1}),

    html.Div(children=[
        html.Label('Checkboxes'),
        dcc.Checklist(all_filenames,),

        html.Br(),
        html.Label('Text Input'),
        dcc.Input(value='MTL', type='text'),

        html.Br(),
        html.Label('Slider'),
        dcc.Slider(
            min=0,
            max=9,
            marks={i: f'Label {i}' if i == 1 else str(i) for i in range(1, 6)},
            value=5,
        )
    ], style={'padding': 10, 'flex': 1})])


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

        
    
