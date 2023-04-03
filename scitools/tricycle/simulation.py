from typing import Any

import numpy as np
import numpy.typing as npt

from scitools.constants import τ


def gen_ellipse(
    locs: npt.NDArray[Any],
    amps: npt.NDArray[Any],
    offset: float | npt.NDArray[np.float32 | np.float64],
    n_time: int = 5000,
) -> npt.NDArray[Any]:
    """
    Generate λ for the negative binomial
    """
    if len(locs) != len(amps):
        raise ValueError("locs and amps must have the same length")

    timepoints = np.linspace(0, τ, n_time)

    λ = amps[:, None] * np.cos(timepoints[None, :] - locs[:, None]) + offset
    return λ.reshape(len(locs), -1)  # (n_genes, n_time)


def simul_negbin(
    locs: npt.NDArray[Any],
    amps: npt.NDArray[Any],
    offset: float,
    bcv: float = 0.1,
    libsize: int = 2000,
    seed: int = 42,
    n_time: int = 5000,
) -> npt.NDArray[np.int_]:
    """
    Simulate negative binomial data
    """
    if len(locs) != len(amps):
        raise ValueError("locs and amps must be the same length")

    n_genes = len(locs)

    λ = gen_ellipse(locs, amps, offset)
    λ /= np.sum(λ, axis=0) / libsize

    rand = np.random.RandomState(seed)
    λ_bcv = rand.gamma(shape=1 / bcv**2, scale=λ * bcv**2, size=(n_genes, n_time))
    counts = rand.poisson(λ_bcv)
    return counts
