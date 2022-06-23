# import functions and text
from matplotlib.pyplot import title
from bin.dashboard.dashboard_functions import read_compsite_files, read_DFI_csvs, read_hng_files
from bin.dashboard.text_boxes import intro_md, flex_md, composite_md, hinges_md
import os
from bin.custom_top import  write_custom_topology
from pathlib import Path

# File management/OS
from pathlib import Path, PurePosixPath


#Plotting
# import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

## DASH
from dash import Dash, dcc, html, Input, Output, State
import subprocess

#Custom topology
from bin.custom_top import make_rb_list



## Config

#external_stylesheets = ["./bWLwgP.css"]
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = Dash(__name__, external_stylesheets=external_stylesheets)
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True
server = app.server
all_outputs = []
for child in Path("output/").iterdir():
    if child.is_dir():
        all_outputs.append(str(child))



### DASH LAYOUT
app.layout = html.Div([
    html.Div(children=[
        html.Label('Select Output'),
        dcc.Dropdown(
            options=all_outputs,
            value=all_outputs[0],
            id="output-dropdown",
            multi=False),

        dcc.Markdown(
            children=intro_md
        ),
        dcc.Graph(
            id='coverage-plot',
        ),
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
    ], style={'padding': 10, 'flex': 1}
    ),

    html.Div(children=[
        dcc.Markdown(
            children=composite_md
        ),
        dcc.Graph(
            id='composite-plot',
        ),
        html.Label('Select Output')
        ,
        dcc.Checklist(
            id="customtop-checklist"
        ),
    ], style={'padding': 10, 'flex': 1}
    ),
    
    html.Div(children=[
        dcc.Markdown(
            children=hinges_md
        ),
    ], style={'padding': 10, 'flex': 1}
    ),
    html.Div(
        id='hinges-output'
    ),

    html.Div(
        id='custom-top-output'
    ),
    html.Div(children=[
        dcc.Input(id='custom-hinges-input', 
            placeholder='Enter hinges here:', 
            value="0:0", type='text') 
    ], style={'padding': 10, 'flex': 1}
    ),

    html.Div(children=[
        html.Button(
            "Generate IMP Topology File",
            id="create-topology-button",
            n_clicks=0),
    ], style={'padding': 10, 'flex': 1}
    ),

    
])


# Update coverage plot
@app.callback(
    Output('coverage-plot', 'figure'),
    [Input('output-dropdown', 'value')]
    )
def update_graph(options_chosen):
    i = 0
    df_list = []
    structure_list = []
    for child in Path(os.path.join(options_chosen, "REPORT", "COVERAGE")).iterdir():
        if child.is_file() and "composite" not in str(child):
            i += 1
            df = pd.read_csv(child)
            df_list.append(df)
            structure_list.append(child)
                
    fig1 = make_subplots(rows=i, cols=1, shared_xaxes=True, x_title="Residue position")

    i = 1
    for df in df_list:
        fig1.append_trace(go.Scatter(
            x=df[df.columns[0]], # ResID
            y=df[df.columns[1]],
            fill='tozeroy', 
            name=str(structure_list[i-1])
        ), row=i, col=1)
        i +=1

    # fig1.update_layout(height=400, width=1000, title_text="Coverage")
    fig1.update_layout(title_text="Coverage")
    fig1.update_yaxes(showgrid=False, range=[0,1], showticklabels=False)    
    return fig1

# Update Flex Plot
@app.callback(
    Output('flex-plot', 'figure'),
    [Input('output-dropdown', 'value')]
    )
def update_graph(options_chosen):
    dfi_dict = read_DFI_csvs(os.path.join(options_chosen, "REPORT", "DFI"))
    hng_dict = read_hng_files(os.path.join(options_chosen, "HINGES"))
    dfi_files = [dfi_file for dfi_file in dfi_dict.keys()  if "AF_DFI" not in str(dfi_file.stem) and "RF_DFI" not in str(dfi_file.stem)]

    fig2 = make_subplots(
        rows=len(dfi_files), cols=1, shared_xaxes=True, 
        x_title="Residue position"
        )

    i = 1
    for dfi_file in dfi_files:
        df = dfi_dict[dfi_file]
        fig2.append_trace(go.Scatter(
            x=df[df.columns[0]], # resIDs
            y=df[df.columns[1]], # pctdfi
            name=str(dfi_file)
        ), row=i, col=1)
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
    fig2.update_layout(title_text="DFI profiles + Predicted hinges", 
                      margin_pad=10, barmode="group", legend=dict(orientation="h",  y=-0.35))
    fig2.update_yaxes(showgrid=False, range=[0,1], nticks=2)
    return fig2



# Update options for checklist for custom topology
@app.callback(
    Output('customtop-checklist', 'options'),
    Output('customtop-checklist', 'value'),
    [Input('output-dropdown', 'value')]
    )
def update_dropdown(selected_output):
    structure_list = []
    for child in Path(os.path.join(selected_output, "REPORT", "COVERAGE")).iterdir():
        if child.is_file() and "composite" not in str(child):
            structure_list.append(str(child))

    
    comp_dict = read_compsite_files(os.path.join(selected_output, "REPORT", "COVERAGE"))
    split_path = selected_output.split("/")
    out_name = split_path[-1]
    filename = os.path.join(selected_output, "REPORT", "COVERAGE", f"{out_name}_composite_coverage.csv")
    df = pd.read_csv(filename)
    initial_structure_list = []
    for str1 in df.columns[1:]:
        id1 = str(os.path.basename(str1))
        for str2 in structure_list:
            id2 = str(os.path.basename(str2))
            print(f"COMPARING: {id1} and {id2}")
            if id1[0:-4] == id2[0:-13]:
                initial_structure_list.append(str2)
                print(f"ADDED {str2}")
    
        
    return structure_list, initial_structure_list


# Update custom composite
@app.callback(
    Output('composite-plot', 'figure'),
    [Input('customtop-checklist', 'value')]
    )
def update_graph(options_chosen):    
    if len(options_chosen) == 0:
        return None

    if options_chosen is None:
        fig = {}
        return fig

    i = 0
    df_list = []
    structure_list = []
    for file in options_chosen:
            if "composite" not in str(file):
                i += 1
                df = pd.read_csv(file)
                df_list.append(df)
                structure_list.append(file)
                
    fig4 = make_subplots(rows=i, cols=1, shared_xaxes=True, x_title="Residue position")

    i = 1
    for df in df_list:
        fig4.append_trace(go.Scatter(
            x=df[df.columns[0]], # ResID
            y=df[df.columns[1]],
            fill='tozeroy', 
            name=str(structure_list[i-1])
        ), row=i, col=1)
        i +=1

    # fig4.update_layout(height=400, width=1000, title_text="Coverage")
    fig4.update_layout(title_text="Coverage")
    fig4.update_yaxes(showgrid=False, range=[0,1], showticklabels=False)

    
    return fig4

# Create topology file on click
@app.callback(
    Output(component_id='hinges-output', component_property='children'),
    State(component_id='customtop-checklist', component_property='value'),
    Input(component_id="create-topology-button", component_property="n_clicks"),
    State(component_id='output-dropdown',  component_property='value'),
    State(component_id="custom-hinges-input", component_property="value")
)
def onclick_topology(selected_fragments, n_clicks, output_dir, str_hinges_input):
    structure_list = []
    clicks = n_clicks
    try:
        for child in Path(os.path.join(output_dir, "PDB", "total")).iterdir():
             if child.is_file() and "composite" not in str(child):
                for name in selected_fragments:
                    if os.path.basename(child)[0:-4] == os.path.basename(name)[0:-13]:
                        structure_list.append(child)
    except:
        pass

    try:
        for child in Path(os.path.join(output_dir, "PDB", "partial")).iterdir():
            if child.is_file() and "composite" not in str(child):
                for name in selected_fragments:
                    if os.path.basename(child)[0:-4] == os.path.basename(name)[0:-13]:
                        structure_list.append(child)
    except:
        pass
    
    try:
        for child in Path(os.path.join(output_dir, "PDB", "CHAINS")).iterdir():
            if child.is_file() and "composite" not in str(child):
                for name in selected_fragments:
                    if os.path.basename(child)[0:-4] == os.path.basename(name)[0:-13]:
                        structure_list.append(child)
    except:
        pass
    
    try:
        for child in Path(os.path.join(output_dir, "ALPHAFOLD", "DOMAINS")).iterdir():
            if child.is_file() and "confident" not in str(child) and "domains" not in str(child):
                for name in selected_fragments:
                    if str(os.path.basename(child)[0:-4]) == str(os.path.basename(name)[0:-13]):
                        structure_list.append(child)
    except:
        pass
    
    try:
        for child in Path(os.path.join(output_dir, "ROSETTAFOLD", "DOMAINS")).iterdir():
            if child.is_file() and "confident" not in str(child) and "domains" not in str(child):
                for name in selected_fragments:
                    if str(os.path.basename(child)[0:-4]) == str(os.path.basename(name)[0:-13]):
                        structure_list.append(child)
    except:
        pass
     
    fasta = "input_fasta/"+str(os.path.basename(output_dir))+".fasta"

    rigid_bodies = make_rb_list(structure_list, fasta)

    ## incorporate the hinges
    hinges_list = [hinge for hinge in str_hinges_input.split(",")] 
    hinges_list = [tuple(i.split(':')) for i in hinges_list]
    hinges_list = [(int(i[0]), int(i[1])) for i in hinges_list]


    final_rigid_bodies = []
    
    for rb in rigid_bodies:
        print(f"RB: {rb.pdb_fn}")
        split_rb = rb.split_rb_hinges(hinges_list)
        print(f"RB SPLIT: {[rb.residue_range for rb in split_rb]}")
        final_rigid_bodies = final_rigid_bodies + split_rb

    print(f"INITIAL = {len(rigid_bodies)}, FINAL = {len(final_rigid_bodies)}")

    
    final_rigid_bodies.sort(key=lambda x: x.residue_range[0])
    str_out = str(output_dir)
    out_name = str_out.split("/")[-1]
    # Write the topology file
    write_custom_topology(os.path.join(output_dir, "IMP", f"{out_name}_custom.topology"), final_rigid_bodies)
    
     
    
    return f"Topology file created with:{[str(rb.pdb_fn) for rb in final_rigid_bodies]}"


if __name__ == "__main__":
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

        
    
