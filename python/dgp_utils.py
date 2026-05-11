from __future__ import annotations
from typing import Dict, Optional, Tuple
import numpy as np
Array = np.ndarray
# ============================================================
# Missingness helpers
# ============================================================

def generate_W_mcar(Y_train: Array, missing_rate: float = 0.30, seed: Optional[int] = None) -> Array:
    """
    MCAR missingness mask.

    Parameters
    ----------
    Y_train : ndarray, shape (N, T)
        Training panel.
    missing_rate : float
        Probability an entry is missing.
    seed : int or None
        Random seed.

    Returns
    -------
    W : ndarray, shape (N, T)
        Binary mask with W[i, t] = 1 if observed and 0 if missing.
    """
    if not (0.0 <= missing_rate < 1.0):
        raise ValueError("missing_rate must lie in [0, 1).")

    rng = np.random.default_rng(seed)
    N, T = Y_train.shape
    W = (rng.uniform(size=(N, T)) > missing_rate).astype(np.float32)
    return W


def apply_W(Y: Array, W: Array) -> Array:
    """Apply a binary observation mask to a panel matrix."""
    if Y.shape != W.shape:
        raise ValueError(f"Y.shape={Y.shape} and W.shape={W.shape} must match.")
    return Y * W


# ============================================================
# Small linear-algebra helpers
# ============================================================

def _spectral_radius(M: Array) -> float:
    return float(np.max(np.abs(np.linalg.eigvals(M))))


def _check_square_matrix(M: Array, name: str) -> None:
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError(f"{name} must be square; got shape {M.shape}.")


# ============================================================
# DGP1: linear factor model with stable VARMA(1,1) factors
# ============================================================

def _simulate_varma11_factors(
    T_total: int,
    A1: Array,
    B1: Array,
    shock_cov: Optional[Array] = None,
    burnin: int = 200,
    seed: Optional[int] = None,
) -> Tuple[Array, Dict[str, Array]]:
    """
    Simulate latent factors from

        z_t = A1 z_{t-1} + u_t + B1 u_{t-1},

    where u_t ~ N(0, shock_cov).

    Returns
    -------
    z : ndarray, shape (T_total, r)
        Latent path after burn-in.
    state : dict
        Contains the latent path, innovation path, and model parameters.
    """
    _check_square_matrix(A1, "A1")
    _check_square_matrix(B1, "B1")

    if A1.shape != B1.shape:
        raise ValueError("A1 and B1 must have the same shape.")

    r = A1.shape[0]

    if _spectral_radius(A1) >= 1.0:
        raise ValueError("A1 is not stable: spectral radius must be < 1.")

    if shock_cov is None:
        shock_cov = 0.04 * np.eye(r)
    if shock_cov.shape != (r, r):
        raise ValueError(f"shock_cov must have shape {(r, r)}.")

    rng = np.random.default_rng(seed)
    T_sim = T_total + burnin

    z = np.zeros((T_sim, r), dtype=float)
    u = rng.multivariate_normal(mean=np.zeros(r), cov=shock_cov, size=T_sim)

    for t in range(1, T_sim):
        z[t] = A1 @ z[t - 1] + u[t] + B1 @ u[t - 1]

    z_out = z[burnin:].copy()
    u_out = u[burnin:].copy()

    state = {
        "z": z_out,
        "u": u_out,
        "A1": A1,
        "B1": B1,
        "shock_cov": shock_cov,
    }
    return z_out, state


def generate_dgp1_linear_varma(
    T_total: int = 515,
    N: int = 64,
    r: int = 2,
    noise_std: float = 0.10,
    seed: int = 123,
) -> Tuple[Array, Array, Dict]:
    """
    DGP1: linear observation model with stable VARMA(1,1) latent dynamics.

        Y_t = Lambda z_t + eps_t
        z_t = A1 z_{t-1} + u_t + B1 u_{t-1}

    Returns
    -------
    Y : ndarray, shape (N, T_total)
    z : ndarray, shape (T_total, r)
    params : dict
    """
    if r != 2:
        raise ValueError("This rebuttal version is currently set up for r=2.")

    rng = np.random.default_rng(seed)

    A1 = np.array([
        [0.58, 0.10],
        [0.00, 0.48],
    ])
    B1 = np.array([
        [0.22, 0.00],
        [0.06, 0.18],
    ])

    z, latent_state = _simulate_varma11_factors(
        T_total=T_total,
        A1=A1,
        B1=B1,
        shock_cov=0.04 * np.eye(r),
        burnin=200,
        seed=seed + 11,
    )

    Lambda = rng.normal(loc=0.0, scale=1.0, size=(N, r))
    eps = rng.normal(loc=0.0, scale=noise_std, size=(N, T_total))
    Y = Lambda @ z.T + eps

    params = {
        "name": "dgp1",
        "type": "linear",
        "obs_type": "linear",
        "latent_type": "varma11",
        "N": N,
        "r": r,
        "Lambda": Lambda,
        "noise_std": noise_std,
        "A1": latent_state["A1"],
        "B1": latent_state["B1"],
        "shock_cov": latent_state["shock_cov"],
        "u_path": latent_state["u"],
    }
    return Y, z, params


# ============================================================
# DGP2 and DGP3 from the earlier notebook naming convention
# ============================================================

def generate_dgp2(
    T_total: int = 515,
    N: int = 64,
    r: int = 2,
    noise_std: float = 0.10,
    seed: int = 123,
) -> Tuple[Array, Array, Dict]:
    """
    DGP2: nonlinear latent dynamics with linear observation model.

    This matches the role of the earlier notebook's DGP2.
    """
    if r != 2:
        raise ValueError("generate_dgp2 is currently set up for r=2.")

    rng = np.random.default_rng(seed)
    z = np.zeros((T_total, r), dtype=float)

    for t in range(1, T_total):
        z_prev = z[t - 1].copy()
        z[t, 0] = 0.60 * z_prev[0] + 0.20 * np.sin(z_prev[1]) + rng.normal(scale=0.15)
        z[t, 1] = 0.50 * z_prev[1] + 0.15 * (z_prev[0] ** 2) + rng.normal(scale=0.15)

    Lambda = rng.normal(size=(N, r))
    eps = rng.normal(scale=noise_std, size=(N, T_total))
    Y = Lambda @ z.T + eps

    params = {
        "name": "dgp2",
        "type": "nonlinear",
        "obs_type": "linear",
        "latent_type": "nonlinear",
        "N": N,
        "r": r,
        "Lambda": Lambda,
        "noise_std": noise_std,
    }
    return Y, z, params


def generate_dgp3(
    T_total: int = 515,
    N: int = 64,
    r: int = 2,
    noise_std: float = 0.10,
    seed: int = 123,
) -> Tuple[Array, Array, Dict]:
    """
    DGP3: nonlinear latent dynamics with nonlinear observation model.

    This matches the role of the earlier notebook's DGP3.
    """
    if r != 2:
        raise ValueError("generate_dgp3 is currently set up for r=2.")

    rng = np.random.default_rng(seed)
    z = np.zeros((T_total, r), dtype=float)

    for t in range(1, T_total):
        z_prev = z[t - 1].copy()
        z[t, 0] = 0.65 * z_prev[0] + 0.15 * np.sin(z_prev[1]) + rng.normal(scale=0.15)
        z[t, 1] = 0.55 * z_prev[1] + 0.10 * (z_prev[0] ** 2) + rng.normal(scale=0.15)

    Wlin = rng.normal(size=(N, r))
    Vnon = rng.normal(size=(N, r))

    Y = np.zeros((N, T_total), dtype=float)
    for t in range(T_total):
        Y[:, t] = (
            Wlin @ z[t]
            + 0.40 * np.sin(Vnon @ z[t])
            + rng.normal(scale=noise_std, size=N)
        )

    params = {
        "name": "dgp3",
        "type": "nonlinear",
        "obs_type": "nonlinear",
        "latent_type": "nonlinear",
        "N": N,
        "r": r,
        "Wlin": Wlin,
        "Vnon": Vnon,
        "noise_std": noise_std,
    }
    return Y, z, params


# ============================================================
# Unified dispatcher
# ============================================================

def generate_panel(
    dgp_name: str,
    T_total: int,
    N: int,
    r: int = 2,
    noise_std: float = 0.10,
    seed: int = 123,
) -> Tuple[Array, Array, Dict]:
    """
    Unified DGP dispatcher.

    Parameters
    ----------
    dgp_name : str
        One of {'dgp1', 'dgp2', 'dgp3'}.
    """
    dgp_name = dgp_name.lower()

    if dgp_name == "dgp1":
        return generate_dgp1_linear_varma(T_total=T_total, N=N, r=r, noise_std=noise_std, seed=seed)
    if dgp_name == "dgp2":
        return generate_dgp2(T_total=T_total, N=N, r=r, noise_std=noise_std, seed=seed)
    if dgp_name == "dgp3":
        return generate_dgp3(T_total=T_total, N=N, r=r, noise_std=noise_std, seed=seed)

    raise ValueError("Unknown dgp_name. Expected one of {'dgp1', 'dgp2', 'dgp3'}.")


# ============================================================
# Conditional mean target C0 = E[Y_{T+1:T+H} | present state]
# ============================================================

def obs_mean_from_z(z: Array, params: Dict) -> Array:
    """Return E[Y_t | z_t] with output shape (N,)."""
    if params["name"] in {"dgp1", "dgp2"}:
        return params["Lambda"] @ z

    if params["name"] == "dgp3":
        return params["Wlin"] @ z + 0.40 * np.sin(params["Vnon"] @ z)

    raise ValueError(f"Unknown DGP name: {params['name']}")


def latent_step_dgp2(z_prev: Array, rng: np.random.Generator) -> Array:
    z_next = np.zeros_like(z_prev)
    z_next[0] = 0.60 * z_prev[0] + 0.20 * np.sin(z_prev[1]) + rng.normal(scale=0.15)
    z_next[1] = 0.50 * z_prev[1] + 0.15 * (z_prev[0] ** 2) + rng.normal(scale=0.15)
    return z_next


def latent_step_dgp3(z_prev: Array, rng: np.random.Generator) -> Array:
    z_next = np.zeros_like(z_prev)
    z_next[0] = 0.65 * z_prev[0] + 0.15 * np.sin(z_prev[1]) + rng.normal(scale=0.15)
    z_next[1] = 0.55 * z_prev[1] + 0.10 * (z_prev[0] ** 2) + rng.normal(scale=0.15)
    return z_next


def analytic_conditional_mean_C0_dgp1(z_train: Array, params: Dict, H: int = 3) -> Array:
    """
    Exact conditional mean target for DGP1.

    Since
        z_{t+1} = A1 z_t + u_{t+1} + B1 u_t,
    the conditional mean given the present state (z_T, u_T) uses future
    innovations with zero mean.
    """
    if params["name"] != "dgp1":
        raise ValueError("analytic_conditional_mean_C0_dgp1 is only for dgp1.")

    z_t = z_train[-1].copy()
    u_t = params["u_path"][z_train.shape[0] - 1].copy()
    A1 = params["A1"]
    B1 = params["B1"]

    latent_means = np.zeros((H, z_t.shape[0]), dtype=float)

    mu_curr = z_t.copy()
    u_curr = u_t.copy()
    for h in range(H):
        mu_next = A1 @ mu_curr + B1 @ u_curr
        latent_means[h] = mu_next
        mu_curr = mu_next
        u_curr = np.zeros_like(u_curr)

    C0 = np.column_stack([obs_mean_from_z(latent_means[h], params) for h in range(H)])
    return C0


def mc_conditional_mean_C0_nonlinear(
    z_T: Array,
    params: Dict,
    H: int = 3,
    n_mc: int = 2000,
    seed: int = 999,
) -> Array:
    """
    Monte Carlo approximation to C0 for DGP2 and DGP3.

    Returns
    -------
    C0 : ndarray, shape (N, H)
    """
    if params["name"] not in {"dgp2", "dgp3"}:
        raise ValueError("mc_conditional_mean_C0_nonlinear is only for dgp2 and dgp3.")

    rng = np.random.default_rng(seed)
    N = params["N"]
    C0 = np.zeros((N, H), dtype=float)

    for _ in range(n_mc):
        z_curr = z_T.copy()
        for h in range(H):
            if params["name"] == "dgp2":
                z_curr = latent_step_dgp2(z_curr, rng)
            else:
                z_curr = latent_step_dgp3(z_curr, rng)
            C0[:, h] += obs_mean_from_z(z_curr, params)

    C0 /= n_mc
    return C0


def get_conditional_mean_C0(
    z_train: Array,
    params: Dict,
    H: int = 3,
    n_mc: int = 2000,
    seed: int = 999,
) -> Array:
    """
    Unified wrapper for the oracle target

        C0 = E[Y_{T+1:T+H} | present state].

    For dgp1 this is analytic.
    For dgp2 and dgp3 this is Monte Carlo.
    """
    if params["name"] == "dgp1":
        return analytic_conditional_mean_C0_dgp1(z_train=z_train, params=params, H=H)

    if params["name"] in {"dgp2", "dgp3"}:
        return mc_conditional_mean_C0_nonlinear(
            z_T=z_train[-1].copy(),
            params=params,
            H=H,
            n_mc=n_mc,
            seed=seed,
        )

    raise ValueError(f"Unknown params['name']: {params['name']}")


# ============================================================
# Tiny smoke test
# ============================================================

def _smoke_test() -> None:
    for dgp_name in ["dgp1", "dgp2", "dgp3"]:
        Y, z, params = generate_panel(dgp_name=dgp_name, T_total=40, N=8, r=2, seed=123)
        assert Y.shape == (8, 40)
        assert z.shape == (40, 2)

        C0 = get_conditional_mean_C0(z_train=z[:30], params=params, H=3, n_mc=200)
        assert C0.shape == (8, 3)

        W = generate_W_mcar(Y[:, :30], missing_rate=0.25, seed=456)
        Y_obs = apply_W(Y[:, :30], W)
        assert W.shape == (8, 30)
        assert Y_obs.shape == (8, 30)

    print("dgp_utils.py smoke test passed.")


if __name__ == "__main__":
    _smoke_test()
