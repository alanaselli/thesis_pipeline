from dash import Dash, html, dcc
import dash_bootstrap_components as dbc
from components import plot_ROH_profile

def tab1(app: Dash) -> dbc.Tab:
    return dbc.Tab(
        label="Tab one",
        children=[
            # first row
            #dbc.Row([html.H1(app.title), html.Hr()]),
            dbc.Row([html.Hr()]),
            # second row
            dbc.Row(
                [
                    # first col
                    dbc.Col([plot_ROH_profile.render(app)])
                ]
            )
        ]
    )

def create_layout(app: Dash) -> dbc.Container:
    return dbc.Container(tab1(app))