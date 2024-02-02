from dash import Dash
from dash_bootstrap_components.themes import BOOTSTRAP, LITERA

from components.layout import create_layout

def main() -> None:
    app = Dash(external_stylesheets=[LITERA])
    app.title = "TITLE"
    app.layout = create_layout(app)
    app.run()

if __name__ == "__main__":
    main()