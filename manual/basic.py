"""
STEP 1 — 기본 ideogram
- 염색체 클릭 → rotated 상태 표시
- 염색체 선택 dropdown
- 이것만 먼저 확인하고 Step 2로 진행

pip install dash dash-bio
python step1_basic.py → http://localhost:8050
"""

from dash import Dash, html, dcc, Input, Output, callback
import dash_bio as dashbio

app = Dash(__name__)

app.layout = html.Div([

    html.H3("Step 1 — 기본 Ideogram", style={"padding": "16px 24px 0"}),

    # 염색체 선택
    html.Div([
        html.Label("표시할 염색체:", style={"fontWeight": "600", "marginBottom": "6px"}),
        dcc.Dropdown(
            id="chrom-select",
            options=[{"label": str(i), "value": str(i)} for i in range(1, 23)]
                  + [{"label": "X", "value": "X"}, {"label": "Y", "value": "Y"}],
            multi=True,
            value=[str(i) for i in range(1, 23)] + ["X"],
        ),
    ], style={"padding": "12px 24px"}),

    # Ideogram — 기본 prop만 사용
    dashbio.Ideogram(
        id="ideo",
        organism="human",
        assembly="GRCh38",
        orientation="vertical",     # 'vertical' | 'horizontal'
        chrHeight=300,
        chrWidth=12,
        chrMargin=10,
        rows=1,
        rotatable=True,             # 클릭 시 해당 염색체 확대
        showBandLabels=True,        # q21, p13 등 band label
        showChromosomeLabels=True,
        resolution=550,             # 450 | 550 | 850
    ),

    # rotated 상태 출력 — 어떤 값이 오는지 확인용
    html.Div(id="rotated-out", style={
        "padding": "12px 24px",
        "fontFamily": "monospace",
        "background": "#f0f4f8",
        "margin": "12px 24px",
        "borderRadius": "6px",
    }),

], style={"fontFamily": "sans-serif", "maxWidth": "1200px"})


@callback(
    Output("ideo", "chromosomes"),
    Input("chrom-select", "value"),
)
def update_chroms(val):
    return val or [str(i) for i in range(1, 23)]


@callback(
    Output("rotated-out", "children"),
    Input("ideo", "rotated"),
)
def show_rotated(rotated):
    return f"rotated = {rotated!r}  ← 염색체 클릭 시 True/False 토글"


if __name__ == "__main__":
    app.run(debug=True, port=8050)