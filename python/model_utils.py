from __future__ import annotations
from typing import Dict, Optional, Tuple
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from statsmodels.tsa.api import VAR
Array = np.ndarray

# ============================================================
# Generic helpers
# ============================================================

def get_default_device(device: Optional[torch.device] = None) -> torch.device:
    if device is not None:
        return device
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _check_panel_shapes(Y_obs: Array, W: Array) -> None:
    if Y_obs.shape != W.shape:
        raise ValueError(f"Y_obs.shape={Y_obs.shape} and W.shape={W.shape} must match.")
    if Y_obs.ndim != 2:
        raise ValueError("Y_obs must be a 2D array of shape (N, T).")


# ============================================================
# PCA with missingness
# Translation of the user's R function into Python
# ============================================================

def est_pca_missing(Y: Array, W: Array, r: int, ridge: float = 1e-6) -> Dict[str, Array]:
    """
    PCA with missingness based on pairwise covariance estimation.

    Parameters
    ----------
    Y : ndarray, shape (N, T)
        Panel with zero-filled missing entries allowed.
    W : ndarray, shape (N, T)
        Binary observation mask. W[i, t] = 1 if observed.
    r : int
        Number of latent factors.
    ridge : float
        Small ridge penalty for unit-wise weighted least squares when needed.

    Returns
    -------
    out : dict
        Contains
        - factor_est : ndarray, shape (T, r)
        - loading_est : ndarray, shape (N, r)
        - Sigma_hat : ndarray, shape (T, T)
        - eigvals : ndarray, shape (r,)
    """
    _check_panel_shapes(Y, W)

    N, T = Y.shape
    if not (1 <= r <= T):
        raise ValueError("r must be between 1 and T.")

    Sigma_hat = np.zeros((T, T), dtype=float)

    for s in range(T):
        for t in range(T):
            idx = np.where((W[:, s] == 1) & (W[:, t] == 1))[0]
            if len(idx) == 0:
                Sigma_hat[s, t] = 0.0
            else:
                Sigma_hat[s, t] = float(np.mean(Y[idx, s] * Y[idx, t]))

    Sigma_hat = 0.5 * (Sigma_hat + Sigma_hat.T)

    eigvals, eigvecs = np.linalg.eigh(Sigma_hat / T)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order][:r]
    eigvecs = eigvecs[:, order][:, :r]

    F_est = np.sqrt(T) * eigvecs  # shape (T, r)
    L_est = np.zeros((N, r), dtype=float)

    for i in range(N):
        obs_idx = np.where(W[i, :] == 1)[0]
        if len(obs_idx) == 0:
            continue

        Fi = F_est[obs_idx, :]                      # n_i x r
        yi = Y[i, obs_idx]                         # n_i
        WFF = Fi.T @ Fi
        WFY = Fi.T @ yi

        try:
            L_est[i, :] = np.linalg.solve(WFF, WFY)
        except np.linalg.LinAlgError:
            L_est[i, :] = np.linalg.solve(WFF + ridge * np.eye(r), WFY)

    return {
        "factor_est": F_est,
        "loading_est": L_est,
        "Sigma_hat": Sigma_hat,
        "eigvals": eigvals,
    }


def reconstruct_from_pca(loadings: Array, factors: Array) -> Array:
    """Return reconstructed panel of shape (N, T) from loadings and factors."""
    return loadings @ factors.T


def decode_pca_forecasts(loadings: Array, z_forecast: Array) -> Array:
    """
    Decode latent forecasts to panel forecasts using the PCA loading matrix.

    Parameters
    ----------
    loadings : ndarray, shape (N, r)
    z_forecast : ndarray, shape (H, r)

    Returns
    -------
    Y_fore : ndarray, shape (N, H)
    """
    return loadings @ z_forecast.T


# ============================================================
# Masked autoencoder from the notebook
# ============================================================

class MaskedAutoencoder(nn.Module):
    def __init__(self, N: int, r: int, hidden: int = 128):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(2 * N, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, r),
        )
        self.decoder = nn.Sequential(
            nn.Linear(r, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, N),
        )

    def encode(self, y: torch.Tensor, W: torch.Tensor) -> torch.Tensor:
        x = torch.cat([y * W, W], dim=-1)
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def forward(self, y: torch.Tensor, W: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        z = self.encode(y, W)
        yhat = self.decode(z)
        return yhat, z


def train_masked_autoencoder(
    Y_obs: Array,
    W: Array,
    r: int = 2,
    hidden: int = 128,
    epochs: int = 300,
    lr: float = 1e-3,
    weight_decay: float = 0.0,
    device: Optional[torch.device] = None,
    verbose: bool = False,
) -> Tuple[MaskedAutoencoder, Array]:
    """
    Train masked autoencoder on zero-filled observed panel.

    Returns
    -------
    model : MaskedAutoencoder
    zhat : ndarray, shape (T, r)
        Encoded latent sequence.
    """
    _check_panel_shapes(Y_obs, W)
    device = get_default_device(device)

    N, T = Y_obs.shape
    model = MaskedAutoencoder(N=N, r=r, hidden=hidden).to(device)

    y_data = torch.tensor(Y_obs.T, dtype=torch.float32, device=device)  # T x N
    w_data = torch.tensor(W.T, dtype=torch.float32, device=device)      # T x N

    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)

    model.train()
    eps = 1e-8
    for epoch in range(epochs):
        optimizer.zero_grad()
        yhat, z = model(y_data, w_data)
        loss = (w_data * (y_data - yhat) ** 2).sum() / (w_data.sum() + eps)
        loss.backward()
        optimizer.step()

        if verbose and ((epoch + 1) % 50 == 0 or epoch == 0):
            print(f"AE epoch {epoch + 1:4d}/{epochs}, loss={loss.item():.6f}")

    model.eval()
    with torch.no_grad():
        zhat = model.encode(y_data, w_data).cpu().numpy()

    return model, zhat


def decode_ae_forecasts(
    ae_model: MaskedAutoencoder,
    z_forecast: Array,
    device: Optional[torch.device] = None,
) -> Array:
    """
    Decode latent forecasts using the AE decoder.

    Parameters
    ----------
    z_forecast : ndarray, shape (H, r)

    Returns
    -------
    Y_fore : ndarray, shape (N, H)
    """
    device = get_default_device(device)
    z_t = torch.tensor(z_forecast, dtype=torch.float32, device=device)

    ae_model.eval()
    with torch.no_grad():
        yhat = ae_model.decode(z_t).cpu().numpy()  # H x N
    return yhat.T


# ============================================================
# Latent dynamics: supervised sequence data
# ============================================================

def build_supervised_latent_data(z: Array, L: int = 5) -> Tuple[Array, Array]:
    """
    Build supervised latent sequences.

    Returns
    -------
    X : ndarray, shape (T-L, L, r)
    Y : ndarray, shape (T-L, r)
    """
    T, r = z.shape
    if L < 1 or L >= T:
        raise ValueError("Need 1 <= L < len(z) to build supervised latent data.")

    X, Y = [], []
    for t in range(L, T):
        X.append(z[t - L:t])
        Y.append(z[t])
    return np.asarray(X), np.asarray(Y)


# ============================================================
# LSTM forecaster from the notebook
# ============================================================

class LatentLSTM(nn.Module):
    def __init__(self, r: int, hidden: int = 64):
        super().__init__()
        self.lstm = nn.LSTM(input_size=r, hidden_size=hidden, batch_first=True)
        self.fc = nn.Linear(hidden, r)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        out, _ = self.lstm(x)
        return self.fc(out[:, -1, :])


def train_latent_lstm(
    zhat: Array,
    L: int = 5,
    hidden: int = 64,
    epochs: int = 200,
    lr: float = 1e-3,
    weight_decay: float = 0.0,
    device: Optional[torch.device] = None,
    verbose: bool = False,
) -> LatentLSTM:
    device = get_default_device(device)

    X, Y = build_supervised_latent_data(zhat, L=L)
    X_t = torch.tensor(X, dtype=torch.float32, device=device)
    Y_t = torch.tensor(Y, dtype=torch.float32, device=device)

    r = zhat.shape[1]
    model = LatentLSTM(r=r, hidden=hidden).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=weight_decay)
    loss_fn = nn.MSELoss()

    model.train()
    for epoch in range(epochs):
        optimizer.zero_grad()
        pred = model(X_t)
        loss = loss_fn(pred, Y_t)
        loss.backward()
        optimizer.step()

        if verbose and ((epoch + 1) % 50 == 0 or epoch == 0):
            print(f"LSTM epoch {epoch + 1:4d}/{epochs}, loss={loss.item():.6f}")

    return model


def forecast_latent_lstm(
    model: LatentLSTM,
    zhat: Array,
    H: int = 3,
    L: int = 5,
    device: Optional[torch.device] = None,
) -> Array:
    device = get_default_device(device)

    history = zhat.copy().tolist()
    preds = []

    model.eval()
    with torch.no_grad():
        for _ in range(H):
            x = np.asarray(history[-L:])[None, :, :]
            x_t = torch.tensor(x, dtype=torch.float32, device=device)
            z_next = model(x_t).cpu().numpy()[0]
            preds.append(z_next)
            history.append(z_next)
    return np.asarray(preds)


# ============================================================
# VAR forecaster for latent factors
# ============================================================

def fit_var_model(
    zhat: Array,
    maxlags: int = 5,
    ic: str = "bic",
    trend: str = "c",
) -> Dict:
    """
    Fit a VAR model on the estimated latent factors.

    Parameters
    ----------
    zhat : ndarray, shape (T, r)
    maxlags : int
    ic : {'aic', 'bic', 'hqic', 'fpe'}
    trend : str

    Returns
    -------
    out : dict
        Contains the fitted statsmodels VARResults and selected lag order.
    """
    if zhat.ndim != 2:
        raise ValueError("zhat must be a 2D array of shape (T, r).")

    T = zhat.shape[0]
    if T <= 3:
        raise ValueError("Need more time points to fit a VAR.")

    maxlags_eff = min(maxlags, max(1, (T - 1) // 3))
    model = VAR(zhat)

    try:
        sel = model.select_order(maxlags=maxlags_eff)
        selected = getattr(sel, ic)
        if selected is None or not np.isfinite(selected) or selected < 1:
            selected = 1
    except Exception:
        selected = 1

    results = model.fit(maxlags=int(selected), ic=None, trend=trend)
    return {
        "var_results": results,
        "selected_lag": int(selected),
        "ic": ic,
    }


def forecast_latent_var(var_results, zhat: Array, H: int = 3) -> Array:
    """
    Forecast latent factors using a fitted VAR model.

    Returns
    -------
    z_fore : ndarray, shape (H, r)
    """
    p = int(var_results.k_ar)
    if zhat.shape[0] < p:
        raise ValueError("zhat is too short for the fitted VAR order.")

    history = zhat[-p:]
    z_fore = var_results.forecast(y=history, steps=H)
    return np.asarray(z_fore)


# ============================================================
# Method wrappers for the ablation study
# ============================================================

def fit_and_forecast_pca_var(
    Y_obs: Array,
    W: Array,
    r: int,
    H: int,
    var_maxlags: int = 5,
    var_ic: str = "bic",
) -> Dict:
    pca_out = est_pca_missing(Y_obs, W, r=r)
    zhat = pca_out["factor_est"]
    loadings = pca_out["loading_est"]

    var_out = fit_var_model(zhat, maxlags=var_maxlags, ic=var_ic)
    z_fore = forecast_latent_var(var_out["var_results"], zhat, H=H)
    Y_fore = decode_pca_forecasts(loadings, z_fore)

    return {
        "method": "pca_var",
        "factor_method": "pca_missing",
        "dyn_method": "var",
        "Y_fore": Y_fore,
        "z_hat": zhat,
        "z_fore": z_fore,
        "loadings": loadings,
        "factor_model": pca_out,
        "dyn_model": var_out,
        "selected_lag": var_out["selected_lag"],
    }


def fit_and_forecast_pca_lstm(
    Y_obs: Array,
    W: Array,
    r: int,
    H: int,
    lookback_L: int = 5,
    lstm_hidden: int = 64,
    lstm_epochs: int = 200,
    lr: float = 1e-3,
    weight_decay: float = 0.0,
    device: Optional[torch.device] = None,
    verbose: bool = False,
) -> Dict:
    pca_out = est_pca_missing(Y_obs, W, r=r)
    zhat = pca_out["factor_est"]
    loadings = pca_out["loading_est"]

    lstm_model = train_latent_lstm(
        zhat,
        L=lookback_L,
        hidden=lstm_hidden,
        epochs=lstm_epochs,
        lr=lr,
        weight_decay=weight_decay,
        device=device,
        verbose=verbose,
    )
    z_fore = forecast_latent_lstm(lstm_model, zhat, H=H, L=lookback_L, device=device)
    Y_fore = decode_pca_forecasts(loadings, z_fore)

    return {
        "method": "pca_lstm",
        "factor_method": "pca_missing",
        "dyn_method": "lstm",
        "Y_fore": Y_fore,
        "z_hat": zhat,
        "z_fore": z_fore,
        "loadings": loadings,
        "factor_model": pca_out,
        "dyn_model": lstm_model,
    }


def fit_and_forecast_ae_var(
    Y_obs: Array,
    W: Array,
    r: int,
    H: int,
    ae_hidden: int = 128,
    ae_epochs: int = 300,
    lr: float = 1e-3,
    weight_decay: float = 0.0,
    var_maxlags: int = 5,
    var_ic: str = "bic",
    device: Optional[torch.device] = None,
    verbose: bool = False,
) -> Dict:
    ae_model, zhat = train_masked_autoencoder(
        Y_obs,
        W,
        r=r,
        hidden=ae_hidden,
        epochs=ae_epochs,
        lr=lr,
        weight_decay=weight_decay,
        device=device,
        verbose=verbose,
    )

    var_out = fit_var_model(zhat, maxlags=var_maxlags, ic=var_ic)
    z_fore = forecast_latent_var(var_out["var_results"], zhat, H=H)
    Y_fore = decode_ae_forecasts(ae_model, z_fore, device=device)

    return {
        "method": "ae_var",
        "factor_method": "ae",
        "dyn_method": "var",
        "Y_fore": Y_fore,
        "z_hat": zhat,
        "z_fore": z_fore,
        "ae_model": ae_model,
        "dyn_model": var_out,
        "selected_lag": var_out["selected_lag"],
    }


def fit_and_forecast_ae_lstm(
    Y_obs: Array,
    W: Array,
    r: int,
    H: int,
    ae_hidden: int = 128,
    lstm_hidden: int = 64,
    ae_epochs: int = 300,
    lstm_epochs: int = 200,
    lookback_L: int = 5,
    lr: float = 1e-3,
    weight_decay: float = 0.0,
    device: Optional[torch.device] = None,
    verbose: bool = False,
) -> Dict:
    ae_model, zhat = train_masked_autoencoder(
        Y_obs,
        W,
        r=r,
        hidden=ae_hidden,
        epochs=ae_epochs,
        lr=lr,
        weight_decay=weight_decay,
        device=device,
        verbose=verbose,
    )

    lstm_model = train_latent_lstm(
        zhat,
        L=lookback_L,
        hidden=lstm_hidden,
        epochs=lstm_epochs,
        lr=lr,
        weight_decay=weight_decay,
        device=device,
        verbose=verbose,
    )
    z_fore = forecast_latent_lstm(lstm_model, zhat, H=H, L=lookback_L, device=device)
    Y_fore = decode_ae_forecasts(ae_model, z_fore, device=device)

    return {
        "method": "ae_lstm",
        "factor_method": "ae",
        "dyn_method": "lstm",
        "Y_fore": Y_fore,
        "z_hat": zhat,
        "z_fore": z_fore,
        "ae_model": ae_model,
        "dyn_model": lstm_model,
    }


def fit_and_forecast_method(
    method_name: str,
    Y_obs: Array,
    W: Array,
    r: int,
    H: int,
    **kwargs,
) -> Dict:
    """
    Unified method dispatcher.

    method_name in {'pca_var', 'pca_lstm', 'ae_var', 'ae_lstm'}
    """
    method_name = method_name.lower()

    if method_name == "pca_var":
        return fit_and_forecast_pca_var(Y_obs=Y_obs, W=W, r=r, H=H, **kwargs)
    if method_name == "pca_lstm":
        return fit_and_forecast_pca_lstm(Y_obs=Y_obs, W=W, r=r, H=H, **kwargs)
    if method_name == "ae_var":
        return fit_and_forecast_ae_var(Y_obs=Y_obs, W=W, r=r, H=H, **kwargs)
    if method_name == "ae_lstm":
        return fit_and_forecast_ae_lstm(Y_obs=Y_obs, W=W, r=r, H=H, **kwargs)

    raise ValueError("Unknown method_name. Expected one of {'pca_var', 'pca_lstm', 'ae_var', 'ae_lstm'}." )


# ============================================================
# Tiny smoke test
# ============================================================

def _smoke_test() -> None:
    rng = np.random.default_rng(123)
    N, T, r, H = 12, 30, 2, 3

    Y = rng.normal(size=(N, T))
    W = (rng.uniform(size=(N, T)) > 0.2).astype(float)
    Y_obs = Y * W

    out1 = fit_and_forecast_pca_var(Y_obs, W, r=r, H=H, var_maxlags=3)
    assert out1["Y_fore"].shape == (N, H)
    assert out1["z_hat"].shape == (T, r)
    assert out1["z_fore"].shape == (H, r)

    out2 = fit_and_forecast_pca_lstm(Y_obs, W, r=r, H=H, lookback_L=4, lstm_epochs=5)
    assert out2["Y_fore"].shape == (N, H)

    out3 = fit_and_forecast_ae_var(Y_obs, W, r=r, H=H, ae_epochs=5, var_maxlags=3)
    assert out3["Y_fore"].shape == (N, H)

    out4 = fit_and_forecast_ae_lstm(Y_obs, W, r=r, H=H, ae_epochs=5, lstm_epochs=5, lookback_L=4)
    assert out4["Y_fore"].shape == (N, H)

    print("model_utils.py smoke test passed.")


if __name__ == "__main__":
    _smoke_test()
