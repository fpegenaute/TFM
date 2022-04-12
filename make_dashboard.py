# import functions and text
from bin.dashboard.dashboard_functions import read_compsite_files, read_DFI_csvs, read_hng_files
from bin.dashboard.text_boxes import intro_md, flex_md, composite_md, custom_md
import os
from bin.custom_top import make_composite

# File management/OS
from pathlib import Path, PurePosixPath


#Plotting
# import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd

## DASH
from dash import Dash, dcc, html, Input, Output, State
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

### DASH LAYOUT
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
    ], style={'padding': 10, 'flex': 1}
    ),
    
    html.Div(children=[
        dcc.Markdown(
            children=custom_md
        ),
        html.Label('Select Output'),
        dcc.Dropdown(
            id="customtop-dropdown",
            multi=True),

        dcc.Graph(
            id='customtop-plot',
        ),
    ], style={'padding': 10, 'flex': 1}
    ),
    
    html.Div(children=[
        html.Button(
            "Generate IMP Topology File",
            id="create-topology",
            n_clicks=0),
    ], style={'padding': 10, 'flex': 1}
    ),

    html.Div(
        id='custom-top-output'
    ),
    
])

from bin.utilities  import get_filename_ext, get_chain_names, get_residue_range
import fnmatch
from bin.custom_top import RigidBody, write_custom_topology
from bin.custom_top import make_composite
from pathlib import Path

@app.callback(
    Output(component_id='custom-top-output', component_property='children'),
    Input(component_id='create-topology', component_property='n_clicks'),
    Input(component_id='customtop-dropdown', component_property='value'),
    Input(component_id='frag-dropdown',  component_property='value')
)
def onclick_topology(nclicks, fragments, output_dir):
    if nclicks > 0:
        i = 1
        rigid_bodies = []
        for csv in fragments: 
            frag_name = os.path.basename(csv) 
            for path in Path(output_dir).rglob('*.pdb'):
                if frag_name[0:5] in  os.path.basename(path): # Provisional, only takes the first letters of the PDB file 
                    filename, extension = get_filename_ext(str(path))
                    
                    # Extract chain name
                    chain_IDs = get_chain_names(path)
                    if len(chain_IDs) > 1 or fnmatch.fnmatch(path, "*AF.pdb"):
                        print(f"""{path}: Assuming a AlphaFold model""")

                        for chain in chain_IDs:
                            # Extract residue range
                            res_range = get_residue_range(str(path), chain=chain)    
                            # Create the RigidBody instance
                            rigid_body = RigidBody(resolution="all",
                            molecule_name= f"{filename}_{chain}", 
                            color="orange" , 
                            fasta_fn="TESTING FASTA FILENAME",
                            # fasta_fn=fasta, 
                            pdb_fn=path, 
                            chain=chain,
                            residue_range=res_range , 
                            rigid_body=i, 
                            super_rigid_body="", 
                            chain_of_super_rigid_bodies="", 
                            bead_size=20,
                            em_residues_per_gaussian=0, 
                            type="AF_model")
                            # Add the rigid body to a list
                            rigid_bodies.append(rigid_body)
                            i +=1
                    elif fnmatch.fnmatch(path, "*RF.pdb"):
                        print(f"""{path}: Assuming a RoseTTaFold model""")

                        for chain in chain_IDs:
                            # Extract residue range
                            res_range = get_residue_range(path, chain=chain)    
                            # Create the RigidBody instance
                            rigid_body = RigidBody(resolution="all",
                            molecule_name= f"{filename}_{chain}", 
                            color="orange" , 
                            fasta_fn="TESTING FASTA FILENAME",
                            # fasta_fn=fasta,  
                            pdb_fn=path, 
                            chain=chain,
                            residue_range=res_range , 
                            rigid_body=i, 
                            super_rigid_body="", 
                            chain_of_super_rigid_bodies="", 
                            bead_size=20,
                            em_residues_per_gaussian=0, 
                            type="RF_model")
                            # Add the rigid body to a list
                            rigid_bodies.append(rigid_body)
                            i +=1
                    elif len(chain_IDs) > 1:
                        print(f"Assuming {path},experimental with multiple chains")
                        for chain in chain_IDs:
                            # Extract residue range
                            res_range = get_residue_range(path, chain=chain)    
                            # Create the RigidBody instance
                            rigid_body = RigidBody(resolution="all",
                            molecule_name= f"{filename}_{chain}", 
                            color="blue" , 
                            fasta_fn="TESTING FASTA FILENAME",
                            # fasta_fn=fasta, 
                            pdb_fn=path, 
                            chain=chain,
                            residue_range=res_range , 
                            rigid_body=i, 
                            super_rigid_body="", 
                            chain_of_super_rigid_bodies="", 
                            bead_size=10,
                            em_residues_per_gaussian=0, 
                            type="experimental")
                            # Add the rigid body to a list
                            rigid_bodies.append(rigid_body)
                            i +=1
                    else:
                        # Extract residue range
                        res_range = get_residue_range(path)    
                        # Create the RigidBody instance
                        rigid_body = RigidBody(resolution="all",
                        molecule_name= filename, 
                        color="blue" , 
                        fasta_fn="TESTING FASTA FILENAME",
                        # fasta_fn=fasta, 
                        pdb_fn=path, 
                        chain=chain_IDs[0],
                        residue_range=res_range , 
                        rigid_body=i, 
                        super_rigid_body="", 
                        chain_of_super_rigid_bodies="", 
                        bead_size=10,
                        em_residues_per_gaussian=0, 
                        type="experimental")

                        # Add the rigid body to a list
                        rigid_bodies.append(rigid_body)
                        i +=1

        ## MAKE THE COMPOSITE
        composite_rb = make_composite(rigid_bodies)


        # Convert to list and sort by the ones who start earlier in the sequence
        composite_rb.sort(key=lambda x: x.residue_range[0])

        # Write the topology file
        write_custom_topology(("button_test.topology"), composite_rb)


    return f"""N clicks: {nclicks}, \nFragments Selected: {fragments}"""


# Update coverage plot
@app.callback(
    Output('coverage-plot', 'figure'),
    [Input('frag-dropdown', 'value')]
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

# Update Flex Plot
@app.callback(
    Output('flex-plot', 'figure'),
    [Input('frag-dropdown', 'value')]
    )
def update_graph(options_chosen):
    dfi_dict = read_DFI_csvs(os.path.join(options_chosen, "REPORT", "DFI"))
    hng_dict = read_hng_files(os.path.join(options_chosen, "HINGES"))

    fig2 = make_subplots(rows=len(dfi_dict.keys()), cols=1, shared_xaxes=True)

    i = 1
    for dfi_file in dfi_dict.keys():
        df = dfi_dict[dfi_file]
        fig2.append_trace(go.Scatter(
            x=df[df.columns[0]], # resIDs
            y=df[df.columns[1]], # pctdfi
            name=str(dfi_file)
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
    return fig2

# Update composite Plot
@app.callback(
    Output('composite-plot', 'figure'),
    [Input('frag-dropdown', 'value')]
    )
def update_graph(options_chosen):
    composite_dir = os.path.join(options_chosen, "REPORT", "COVERAGE")
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

    return fig3

# Update options for dropdown for custom topology
@app.callback(
    Output('customtop-dropdown', 'options'),
    [Input('frag-dropdown', 'value')]
    )
def update_dropdown(selected_output):
    structure_list = []
    for child in Path(os.path.join(selected_output, "REPORT", "COVERAGE")).iterdir():
        if child.is_file() and "composite" not in str(child):
            structure_list.append(str(child))
    print(f"STRUCT_LIST {structure_list}")
    # return [{"label": x, "value": x} for x in structure_list]
    return structure_list

# Update values for dropdown custom topology
@app.callback(
    Output('customtop-dropdown', 'value'),
    [Input('customtop-dropdown', 'options')]
    )
def update_dropdown(selected_options):
    return selected_options


# Update custom composite
@app.callback(
    Output('customtop-plot', 'figure'),
    [Input('customtop-dropdown', 'value')]
    )
def update_graph(options_chosen):    
    if len(options_chosen) == 0:
        return None

    i = 0
    df_list = []
    structure_list = []
    for file in options_chosen:
            if "composite" not in str(file):
                i += 1
                df = pd.read_csv(file)
                df_list.append(df)
                structure_list.append(file)
                
    fig4 = make_subplots(rows=i, cols=1, shared_xaxes=True)

    i = 1
    for df in df_list:
        fig4.append_trace(go.Scatter(
            x=df[df.columns[0]], # ResID
            y=df[df.columns[1]],
            fill='tozeroy', 
            name=str(structure_list[i-1])
        ), row=i, col=1)
        i +=1

    fig4.update_layout(height=400, width=1000, title_text="Coverage")
    fig4.update_yaxes(showgrid=False, range=[0,1], nticks=2)
    
    return fig4

# @app.callback(
#     Output('container-button-basic', 'children'),
#     Input('submit-val', 'n_clicks'),
#     State('input-on-submit', 'value')
# )
# def update_output(n_clicks, value):
#     return 'The input value was "{}" and the button has been clicked {} times'.format(
#         value,
#         n_clicks
#     )


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

        
    
