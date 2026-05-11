from __future__ import annotations
from typing import Dict, List, Optional
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
Array = np.ndarray

# ============================================================
# Core error helpers
# ============================================================

def squared_error_matrix(Y_fore: Array, C0: Array) -> Array:
    """
    Elementwise squared forecast error matrix.

    Parameters
    ----------
    Y_fore : ndarray, shape (N, H)
        Forecasted panel values.
    C0 : ndarray, shape (N, H)
        Oracle conditional mean target.

    Returns
    -------
    sqerr : ndarray, shape (N, H)
        Squared error matrix.
    """
    if Y_fore.shape != C0.shape:
        raise ValueError(f"Y_fore.shape={Y_fore.shape} and C0.shape={C0.shape} must match.")
    return (Y_fore - C0) ** 2


def msfe_first_k_units(Y_fore: Array, C0: Array, k: int = 32) -> float:
    """
    Mean squared forecast error over the first k units and all horizons.
    """
    if Y_fore.shape != C0.shape:
        raise ValueError(f"Y_fore.shape={Y_fore.shape} and C0.shape={C0.shape} must match.")

    k_eff = min(k, Y_fore.shape[0])
    sqerr = squared_error_matrix(Y_fore[:k_eff, :], C0[:k_eff, :])
    return float(np.mean(sqerr))


def sd_squared_error_first_k_units(Y_fore: Array, C0: Array, k: int = 32, ddof: int = 1) -> float:
    """
    Standard deviation of squared forecast errors over the first k units
    and all horizons.
    """
    if Y_fore.shape != C0.shape:
        raise ValueError(f"Y_fore.shape={Y_fore.shape} and C0.shape={C0.shape} must match.")

    k_eff = min(k, Y_fore.shape[0])
    sqerr = squared_error_matrix(Y_fore[:k_eff, :], C0[:k_eff, :]).reshape(-1)

    if sqerr.size <= 1:
        return 0.0
    return float(np.std(sqerr, ddof=ddof))


# ============================================================
# Evaluation wrapper for one fitted method
# ============================================================

def evaluate_forecast(
    Y_fore: Array,
    C0: Array,
    k_eval: int = 32,
    ddof: int = 1,
) -> Dict:
    """
    Evaluate one forecast matrix against the oracle conditional mean target.

    Parameters
    ----------
    Y_fore : ndarray, shape (N, H)
    C0 : ndarray, shape (N, H)
    k_eval : int
        Number of units to evaluate from the top of the panel.
    ddof : int
        Degrees of freedom for standard deviation.

    Returns
    -------
    out : dict
        Contains MSFE, SD of squared errors, and raw squared errors.
    """
    if Y_fore.shape != C0.shape:
        raise ValueError(f"Y_fore.shape={Y_fore.shape} and C0.shape={C0.shape} must match.")

    k_eff = min(k_eval, Y_fore.shape[0])

    sqerr_full = squared_error_matrix(Y_fore, C0)
    sqerr_eval = sqerr_full[:k_eff, :]

    return {
        "k_eval": k_eff,
        "msfe_first_k": float(np.mean(sqerr_eval)),
        "sd_sqerr_first_k": float(np.std(sqerr_eval.reshape(-1), ddof=ddof)) if sqerr_eval.size > 1 else 0.0,
        "sqerr_first_k": sqerr_eval,
        "sqerr_full": sqerr_full,
    }


def evaluate_method_output(
    method_out: Dict,
    C0: Array,
    k_eval: int = 32,
    ddof: int = 1,
) -> Dict:
    """
    Evaluate the output dictionary returned by model_utils.fit_and_forecast_*.
    """
    if "Y_fore" not in method_out:
        raise KeyError("method_out must contain key 'Y_fore'.")

    eval_out = evaluate_forecast(
        Y_fore=method_out["Y_fore"],
        C0=C0,
        k_eval=k_eval,
        ddof=ddof,
    )

    combined = dict(method_out)
    combined.update(eval_out)
    return combined

def compare_methods_msfe(df, method_a, method_b):
    """
    Paired comparison between two methods across replicates.
    """
    merged = df[df["method"] == method_a].merge(
        df[df["method"] == method_b],
        on=["dgp", "T_train", "rep"],
        suffixes=("_a", "_b"),
    )

    results = []

    for (dgp, T_train), sub in merged.groupby(["dgp", "T_train"]):
        diff = sub["msfe_first_k_a"] - sub["msfe_first_k_b"]

        mean_diff = diff.mean()
        sd_diff = diff.std(ddof=1)

        try:
            stat, pval = wilcoxon(diff)
        except:
            pval = np.nan

        results.append({
            "dgp": dgp,
            "T_train": T_train,
            "method_a": method_a,
            "method_b": method_b,
            "mean_diff": mean_diff,
            "sd_diff": sd_diff,
            "p_value": pval,
        })

    return pd.DataFrame(results)

# ============================================================
# Row builders for simulation output
# ============================================================

def make_result_row(
    dgp_name: str,
    method_name: str,
    rep: int,
    T_train: int,
    H: int,
    missing_rate: float,
    method_out: Dict,
    eval_out: Dict,
    extra_info: Optional[Dict] = None,
) -> Dict:
    """
    Create one flat result row for later tabulation.
    """
    row = {
        "dgp": dgp_name,
        "method": method_name,
        "rep": int(rep),
        "T_train": int(T_train),
        "H": int(H),
        "missing_rate": float(missing_rate),
        "k_eval": int(eval_out["k_eval"]),
        "msfe_first_k": float(eval_out["msfe_first_k"]),
        "sd_sqerr_first_k": float(eval_out["sd_sqerr_first_k"]),
    }

    if "factor_method" in method_out:
        row["factor_method"] = method_out["factor_method"]
    if "dyn_method" in method_out:
        row["dyn_method"] = method_out["dyn_method"]
    if "selected_lag" in method_out:
        row["selected_lag"] = method_out["selected_lag"]

    if extra_info is not None:
        row.update(extra_info)

    return row


def rows_to_dataframe(rows: List[Dict]) -> pd.DataFrame:
    """
    Convert a list of result rows to a pandas DataFrame.
    """
    if len(rows) == 0:
        return pd.DataFrame()
    return pd.DataFrame(rows)


# ============================================================
# Summary tables
# ============================================================

def summarize_results(df: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize replicate-level results by DGP, method, and training size.

    Output columns:
    - mean_msfe
    - sd_msfe_across_reps
    - mean_sd_sqerr
    - sd_sd_sqerr_across_reps
    - n_reps
    """
    needed = {"dgp", "method", "T_train", "msfe_first_k", "sd_sqerr_first_k"}
    missing = needed.difference(df.columns)
    if missing:
        raise ValueError(f"DataFrame is missing required columns: {sorted(missing)}")

    summary = (
        df.groupby(["dgp", "method", "T_train"], as_index=False)
        .agg(
            mean_msfe=("msfe_first_k", "mean"),
            sd_msfe_across_reps=("msfe_first_k", "std"),
            mean_sd_sqerr=("sd_sqerr_first_k", "mean"),
            sd_sd_sqerr_across_reps=("sd_sqerr_first_k", "std"),
            n_reps=("rep", "count"),
        )
        .sort_values(["dgp", "T_train", "method"])
        .reset_index(drop=True)
    )
    return summary


def summarize_results_with_methods(df: pd.DataFrame) -> pd.DataFrame:
    """
    Same as summarize_results, but retains factor/dynamics method columns if present.
    """
    needed = {"dgp", "method", "T_train", "msfe_first_k", "sd_sqerr_first_k"}
    missing = needed.difference(df.columns)
    if missing:
        raise ValueError(f"DataFrame is missing required columns: {sorted(missing)}")

    group_cols = ["dgp", "method", "T_train"]
    if "factor_method" in df.columns:
        group_cols.append("factor_method")
    if "dyn_method" in df.columns:
        group_cols.append("dyn_method")

    summary = (
        df.groupby(group_cols, as_index=False)
        .agg(
            mean_msfe=("msfe_first_k", "mean"),
            sd_msfe_across_reps=("msfe_first_k", "std"),
            mean_sd_sqerr=("sd_sqerr_first_k", "mean"),
            sd_sd_sqerr_across_reps=("sd_sqerr_first_k", "std"),
            n_reps=("rep", "count"),
        )
        .sort_values(["dgp", "T_train", "method"])
        .reset_index(drop=True)
    )
    return summary


# ============================================================
# Tiny smoke test
# ============================================================

def _smoke_test() -> None:
    rng = np.random.default_rng(123)

    N, H = 10, 3
    C0 = rng.normal(size=(N, H))
    Y_fore = C0 + 0.2 * rng.normal(size=(N, H))

    eval_out = evaluate_forecast(Y_fore, C0, k_eval=6)
    assert "msfe_first_k" in eval_out
    assert "sd_sqerr_first_k" in eval_out
    assert eval_out["sqerr_first_k"].shape == (6, H)

    method_out = {
        "method": "pca_var",
        "factor_method": "pca_missing",
        "dyn_method": "var",
        "Y_fore": Y_fore,
        "selected_lag": 2,
    }
    combined = evaluate_method_output(method_out, C0, k_eval=6)
    assert combined["sqerr_full"].shape == (N, H)

    rows = []
    for rep in range(4):
        rows.append(
            make_result_row(
                dgp_name="dgp1",
                method_name="pca_var",
                rep=rep,
                T_train=32,
                H=H,
                missing_rate=0.3,
                method_out=method_out,
                eval_out=eval_out,
            )
        )

    df = rows_to_dataframe(rows)
    summary = summarize_results(df)

    assert df.shape[0] == 4
    assert summary.shape[0] == 1

    print("eval_utils.py smoke test passed.")


if __name__ == "__main__":
    _smoke_test()