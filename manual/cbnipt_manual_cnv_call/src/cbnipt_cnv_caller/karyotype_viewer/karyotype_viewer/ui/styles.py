"""
Shared CSS-in-Python style dicts and Dash layout helpers.
Import from here rather than duplicating in each Dash app.
"""

from __future__ import annotations
from typing import Any, Optional

from dash import html


# ---------------------------------------------------------------------------
# Base style tokens
# ---------------------------------------------------------------------------
PAGE: dict = {
    "fontFamily": "Inter,Arial,sans-serif",
    "background":  "#EDF2F7",
    "color":       "#2D3748",
    "padding":     "8px",
    "margin":      "0",
    "fontSize":    "12px",
    "boxSizing":   "border-box",
}

CARD: dict = {
    "background":    "white",
    "border":        "1px solid #D9E2EC",
    "borderRadius":  "7px",
    "marginBottom":  "8px",
    "overflow":      "visible",
}

DASHBAR_STYLE: dict = {
    "height":       "30px",
    "boxSizing":    "border-box",
    "padding":      "5px 10px",
    "display":      "flex",
    "alignItems":   "center",
    "gap":          "8px",
    "borderBottom": "1px solid #E2E8F0",
    "background":   "#F7FAFC",
    "fontSize":     "10px",
    "fontWeight":   "700",
    "letterSpacing":".08em",
    "color":        "#718096",
}


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------
def dashbar(title: str, right: Optional[Any] = None) -> html.Div:
    children: list[Any] = [
        html.Span("● ● ●", style={
            "fontSize": "7px", "color": "#CBD5E0", "letterSpacing": "2px"
        }),
        html.Span(title),
    ]
    if right is not None:
        children.append(html.Div(
            right,
            style={"marginLeft": "auto", "display": "flex",
                   "alignItems": "center", "gap": "8px"},
        ))
    return html.Div(children, style=DASHBAR_STYLE)


def card(
    title: str,
    body: Any,
    right: Optional[Any] = None,
    body_style: Optional[dict] = None,
) -> html.Div:
    style: dict = {"padding": "10px 12px", "boxSizing": "border-box"}
    if body_style:
        style.update(body_style)
    return html.Div(
        [dashbar(title, right), html.Div(body, style=style)],
        style=CARD,
    )


# ---------------------------------------------------------------------------
# Micro UI helpers
# ---------------------------------------------------------------------------
def _label(text: str) -> html.Div:
    """Small muted section label."""
    return html.Div(text, style={
        "fontSize": "9px", "fontWeight": "600",
        "letterSpacing": ".06em", "color": "#A0AEC0",
        "textTransform": "uppercase",
    })


def badge(text: str, color: str = "#3182CE") -> html.Span:
    """Monospace pill badge."""
    return html.Span(text, style={
        "fontFamily": "monospace",
        "fontSize": "11px",
        "fontWeight": "700",
        "color": color,
        "background": f"{color}18",
        "border": f"1px solid {color}44",
        "borderRadius": "4px",
        "padding": "2px 7px",
        "display": "inline-block",
    })


def _rgba(hex_color: str, alpha: float = 0.1) -> str:
    """Convert '#RRGGBB' → 'rgba(r,g,b,alpha)'."""
    h = hex_color.lstrip("#")
    if len(h) == 3:
        h = "".join(c * 2 for c in h)
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f"rgba({r},{g},{b},{alpha})"
