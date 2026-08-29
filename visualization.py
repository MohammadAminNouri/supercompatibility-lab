from __future__ import annotations

import numpy as np
import plotly.graph_objects as go

from .core import canonical_plane, cmc_degeneracy


def _plane_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = normal / np.linalg.norm(normal)
    helper = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.85 else np.array([0.0, 1.0, 0.0])
    u = np.cross(n, helper)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    v /= np.linalg.norm(v)
    return u, v


def cmc_surface_figure(cmc: np.ndarray, extent: float = 2.0) -> go.Figure:
    deg = cmc_degeneracy(cmc)
    evals, evecs = np.linalg.eigh(cmc)
    fig = go.Figure()

    if deg.habit_planes:
        grid = np.linspace(-extent, extent, 17)
        uu, vv = np.meshgrid(grid, grid)
        for i, p in enumerate(deg.habit_planes):
            n = canonical_plane(p)
            b1, b2 = _plane_basis(n)
            xyz = uu[..., None] * b1 + vv[..., None] * b2
            fig.add_trace(
                go.Surface(
                    x=xyz[:, :, 0],
                    y=xyz[:, :, 1],
                    z=xyz[:, :, 2],
                    opacity=0.55,
                    showscale=False,
                    name=f"Habit plane {i+1}",
                    hovertemplate="x=%{x:.2f}<br>y=%{y:.2f}<br>z=%{z:.2f}<extra></extra>",
                )
            )
        title = "CMC degeneracy: compatible plane solution"
    else:
        signs = np.sign(evals)
        pos = np.flatnonzero(evals > 1e-9)
        neg = np.flatnonzero(evals < -1e-9)
        if (len(pos), len(neg)) in {(1, 2), (2, 1)}:
            majority = neg if len(neg) == 2 else pos
            single = pos[0] if len(pos) == 1 else neg[0]
            theta = np.linspace(0, 2 * np.pi, 90)
            r = np.linspace(0.0, extent, 45)
            rr, tt = np.meshgrid(r, theta)
            coords_base = np.zeros((theta.size, r.size, 3))
            coords_base[..., majority[0]] = rr * np.cos(tt) / np.sqrt(abs(evals[majority[0]]))
            coords_base[..., majority[1]] = rr * np.sin(tt) / np.sqrt(abs(evals[majority[1]]))
            for sign in (-1.0, 1.0):
                coords = coords_base.copy()
                coords[..., single] = sign * rr / np.sqrt(abs(evals[single]))
                xyz = np.einsum("ij,abj->abi", evecs, coords)
                fig.add_trace(
                    go.Surface(
                        x=xyz[:, :, 0],
                        y=xyz[:, :, 1],
                        z=xyz[:, :, 2],
                        opacity=0.58,
                        showscale=False,
                        name="CMC zero-change cone",
                        hovertemplate="direction ray<br>x=%{x:.2f}<br>y=%{y:.2f}<br>z=%{z:.2f}<extra></extra>",
                    )
                )
            title = "CMC zero-change directions: double cone"
        else:
            fig.add_annotation(text="The zero set is not a real double-cone for these inputs.", x=0.5, y=0.5, showarrow=False)
            title = "CMC zero set"

    fig.add_trace(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers", marker={"size": 4}, name="origin"))
    fig.update_layout(
        title=title,
        margin=dict(l=0, r=0, t=48, b=0),
        height=520,
        scene=dict(
            xaxis_title="Austenite direction x",
            yaxis_title="Austenite direction y",
            zaxis_title="Austenite direction z",
            aspectmode="cube",
        ),
        legend=dict(orientation="h"),
    )
    return fig
