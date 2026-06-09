"""Marginal calibration (protocol Section 9.1, G-SIM check 1).

Fits the simulator's structural marginal-draw parameters to the empirical FAERS
drug- and event-margin histograms from the locked Arm-1 reference frame
(``results/perpair_arm1.parquet`` — the frame that defines the deposited anchor
counts 52 / 44 / 6). The DGP draws FAERS marginals lognormal on the count scale
and derives DAEN marginals by rho-scaling + Poisson resampling, so the fitted
quantities are:

- ``log_mu_drug``  / ``log_sd_drug``  : MLE of the lognormal on per-pair n_drug_F
- ``log_mu_event`` / ``log_sd_event`` : MLE of the lognormal on per-pair n_event_F
- ``N_F``  : the FAERS surviving-case total (a single value across the frame)
- ``rho``  : empirical DAEN/FAERS size ratio (N_D / N_F)
- ``pi``   : empirical true-null fraction (mean of groundTruth == 0)
- ``n_pairs`` : the reference-set size (rows in the anchor frame)

The fit is written to ``results/marginal_calibration.json`` (the deposited
calibration record). ``sim_config.yaml`` is updated separately, by hand, so the
file's integrity-note comments are preserved; this script prints the exact block
to paste. Nothing here is presented as calibrated until that config edit lands
and ``marginals.calibrated`` flips true.

Run: ``python -m simulation.calibrate``  (add ``--frame PATH`` to override).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FRAME = ROOT / "results" / "perpair_arm1.parquet"
OUT_JSON = ROOT / "results" / "marginal_calibration.json"


def _lognormal_mle(counts: np.ndarray) -> tuple[float, float, int]:
    """Lognormal MLE on the count scale: mu, sigma = mean/sd of log(counts).

    Zeros (a marginal of 0 reports) cannot be log-transformed and are dropped;
    the count of dropped zeros is returned so the caller can report coverage.
    """
    x = np.asarray(counts, dtype=float)
    pos = x[x > 0]
    n_zero = int((x <= 0).sum())
    lx = np.log(pos)
    # MLE uses the population SD (ddof=0); the difference vs ddof=1 is negligible
    # at n~500 and the MLE is the deposited convention.
    return float(lx.mean()), float(lx.std(ddof=0)), n_zero


def calibrate(frame: Path = DEFAULT_FRAME) -> dict:
    df = pd.read_parquet(frame)

    mu_drug, sd_drug, z_drug = _lognormal_mle(df["n_drug_F"].to_numpy())
    mu_event, sd_event, z_event = _lognormal_mle(df["n_event_F"].to_numpy())

    N_F = int(df["N_F"].mode().iloc[0])
    N_D = int(df["N_D"].mode().iloc[0])
    rho = N_D / N_F
    pi = float((df["groundTruth"] == 0).mean())
    n_pairs = int(len(df))

    # Anchor counts as they stand in this frame (re-derived, not trusted blindly).
    h1 = df[df["faers_signal"] & df["daen_powered"]
            & df["cross_db_class"].isin(["faers_only", "concordant_positive"])]
    anchor = {
        "daen_powered": int(df["daen_powered"].sum()),
        "h1_stratum": int(len(h1)),
        "faers_only": int((h1["cross_db_class"] == "faers_only").sum()),
    }

    return {
        "calibrated": True,
        "source_frame": str(frame.relative_to(ROOT)),
        "fitted_at": time.strftime("%Y-%m-%d"),
        "n_pairs": n_pairs,
        "marginals": {
            "log_mu_drug": round(mu_drug, 4), "log_sd_drug": round(sd_drug, 4),
            "log_mu_event": round(mu_event, 4), "log_sd_event": round(sd_event, 4),
            "n_zero_drug": z_drug, "n_zero_event": z_event,
        },
        "population": {"N_F": N_F, "N_D": N_D, "rho": round(rho, 5)},
        "pi_empirical": round(pi, 4),
        "anchor_counts_in_frame": anchor,
    }


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description="Fit simulator marginal calibration.")
    p.add_argument("--frame", default=str(DEFAULT_FRAME))
    p.add_argument("--write", action="store_true",
                   help="write results/marginal_calibration.json")
    args = p.parse_args(argv)

    cal = calibrate(Path(args.frame))
    print(json.dumps(cal, indent=2))
    print("\n--- paste into sim_config.yaml (preserving comments) ---")
    m = cal["marginals"]
    print(f"  log_mu_drug:  {m['log_mu_drug']}")
    print(f"  log_sd_drug:  {m['log_sd_drug']}")
    print(f"  log_mu_event: {m['log_mu_event']}")
    print(f"  log_sd_event: {m['log_sd_event']}")
    print(f"  calibrated:   true")
    print(f"population.N_F: {cal['population']['N_F']}   "
          f"(rho={cal['population']['rho']}, pi={cal['pi_empirical']}, "
          f"n_pairs={cal['n_pairs']})")

    if args.write:
        OUT_JSON.write_text(json.dumps(cal, indent=2))
        print(f"\nwrote {OUT_JSON.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
