# %%
from functools import cache
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import scanpy as sc
import seaborn as sns
from matplotlib.axes import Axes
from sklearn.decomposition import PCA

seed = 42
π = np.pi
τ = 2 * π
n_time = 5000
timepoints = np.linspace(0, τ, n_time)
sns.set()


rand = np.random.RandomState(seed)
n_genes = 100


def gen_plot(i: int, axs: list[Axes]):
    arr = gen_simul(i)
    adata = ad.AnnData(arr, dtype=int)
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata, copy=False)
    normed = adata.X

    pca = PCA(n_components=2)
    proj_genes = pca.fit_transform(normed).T
    proj_cells = pca.fit_transform(normed.T).T
    print(pca.explained_variance_ratio_)

    α, s = 0.3, 1
    ax_kwargs = dict(alpha=α, s=s)

    axs[0].scatter(proj_genes[0], proj_genes[1], **ax_kwargs)
    axs[0].set_title("Genes")
    axs[1].scatter(proj_cells[0], proj_cells[1], **ax_kwargs)
    axs[1].set_title("Cell")
    axs[2].scatter(timepoints, proj_cells[0], **ax_kwargs)
    axs[2].set_title("Cell PC1")
    axs[3].scatter(timepoints, proj_cells[1], **ax_kwargs)
    axs[3].set_title("Cell PC2")


fig, axs = plt.subplots(4, 4, figsize=(10, 10))
fig.dpi = 300
libsize = [2000, 500, 100, 10]
for axs, i in zip(axs, libsize):
    gen_plot(i, axs)

plt.tight_layout()
