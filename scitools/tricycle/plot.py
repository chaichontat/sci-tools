from typing import Any, NamedTuple

import cmocean  # colormap, do not remove
import colorcet as cc  # colormap, do not remove
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.axes import Axes
from matplotlib.figure import Figure


def plot_wheel(
    pcs: npt.NDArray[np.float_],
    cmap: str = "cet_colorwheel",
    fig: Figure | None = None,
    colorize_background: bool = True,
    **kwargs: Any,
) -> tuple[Figure, NamedTuple]:

    θ = np.arctan2(pcs[:, 0], pcs[:, 1])
    fig = fig or plt.figure(1, figsize=(6, 6))
    axs = NamedTuple("axs", scatter=Axes, histx=Axes, histy=Axes)(
        scatter=plt.subplot2grid((6, 6), (1, 0), colspan=5, rowspan=5),
        histx=plt.subplot2grid((6, 6), (0, 0), colspan=5),
        histy=plt.subplot2grid((6, 6), (1, 5), rowspan=5),
    )

    H, xedges, yedges = np.histogram2d(pcs[:, 0], pcs[:, 1], bins=(100, 100))
    axs.scatter.axhline(0, color="black", alpha=0.5, markeredgewidth=0, linewidth=0.75)
    axs.scatter.axvline(0, color="black", alpha=0.5, markeredgewidth=0, linewidth=0.75)
    sc = axs.scatter.scatter(
        pcs[:, 0],
        pcs[:, 1],
        alpha=1,
        s=5,
        c=θ,
        cmap=cmap,
        zorder=2,
        linewidth=0.2,
        **kwargs
    )
    sc.set_edgecolor((0.4, 0.4, 0.4, 0.5))
    axs.scatter.set(xlabel="PC1", ylabel="PC2")

    # Background color wheel
    if colorize_background:
        xlim, ylim = axs.scatter.get_xlim(), axs.scatter.get_ylim()
        re, im = np.mgrid[xlim[0] : xlim[1] : 100j, ylim[0] : ylim[1] : 100j]
        n = re + 1j * im
        angle = np.arctan2(n.real, n.imag)
        axs.scatter.pcolormesh(re, im, angle, cmap=cmap, alpha=0.5)

    axs.histx.hist(pcs[:, 0], bins=xedges, alpha=0.8, edgecolor=None)
    axs.histy.hist(pcs[:, 1], bins=yedges, alpha=0.8, orientation="horizontal", edgecolor=None)

    nullfmt = mpl.ticker.NullFormatter()
    axs.histx.xaxis.set_major_formatter(nullfmt)
    axs.histy.yaxis.set_major_formatter(nullfmt)

    fig.tight_layout()
    return fig, axs
