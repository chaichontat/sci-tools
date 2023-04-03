# %%
from functools import cache

import numpy as np
import numpy.typing as npt
import seaborn as sns
from numpy.random import RandomState
from sklearn.decomposition import PCA

from scitools.constants import τ
from scitools.tricycle.plot import plot_wheel
from scitools.tricycle.simulation import simul_negbin

sns.set()


@cache
def gen_simul(rand: RandomState, libsize: int, n_genes: int = 100) -> npt.NDArray[np.int_]:
    return simul_negbin(
        rand.uniform(0, τ, n_genes),
        3 * np.ones(n_genes),
        offset=5,
        bcv=0.1,
        libsize=libsize,
    )


rand = RandomState(42)
sim = gen_simul(rand, 2000)

pca = PCA(n_components=2)
proj_cells = pca.fit_transform(sim.T)

plot_wheel(proj_cells)

# %%
