from dash import Dash, dcc, html
from dash.dependencies import Input, Output

import pandas as pd
import os
import plotly.express as px
import plotly.graph_objects as go

sim_path = os.path.join("../","simulations","50_50")
scenario_01 = os.path.join("scenario_01")

cons_runs = pd.read_csv(os.path.join(sim_path,scenario_01,"sc_01_66_consecutiveRuns.txt"),
                        sep=" ")

candidates_metrics = pd.read_csv(os.path.join(sim_path,scenario_01,"candidates_metrics.txt"),
                        sep=" ")

cons_runs['length_Mps'] = cons_runs['lengthBps']/1000000
cons_runs['from_Mps'] = cons_runs['from']/1000000
cons_runs['to_Mps'] = cons_runs['to']/1000000

def render(app: Dash) -> html.Div:
    select_chrom = 1
    df = cons_runs.loc[cons_runs['chrom']==select_chrom].copy()

    df = df.merge(candidates_metrics.loc[candidates_metrics["gen"]==66,["ID","Fg"]],
                left_on="id",right_on="ID",how="left")
    df['Fg_cat'] = pd.cut(df['Fg'], bins=3, 
                        labels=['Low Fg','Medium Fg','High Fg'])
    df['Fg_cat'] = df['Fg_cat'].astype("str")

    df['id'] = df['id'].astype("str")

    fig = px.timeline(df, x_start="from", x_end="to", y="id", color="Fg_cat",
                    color_discrete_sequence=["red", "blue", "green"])
    fig.update_yaxes(autorange="reversed")

    fig.layout.xaxis.type = 'linear'
    for d in fig.data:
        filt = df['Fg_cat'] == d.name
        d.x = df[filt]['lengthBps'].tolist()

    return html.Div(dcc.Graph(figure=fig), id='PLOT')