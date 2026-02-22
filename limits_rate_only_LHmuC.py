# python limits_rate_only_LHmuC.py --save --outdir plots_limits_LHmuC

import numpy as np
import math
import argparse
from pathlib import Path
import matplotlib.pyplot as plt

# Optional CMS-like style
try:
    import mplhep as hep
    try:
        hep.style.use("CMS")
    except Exception:
        pass
except Exception:
    pass


# ============================================================
# LHmuC inputs (pb) from your MG banner tables
# ============================================================

sigma_SM = 0.5737256299

# CP-even: sigma(TR) points, TI=0
TR_pts = np.array([-0.1, -0.01, 0.0, +0.01, +0.1], dtype=float)
sig_TR = np.array([
    0.4127357714,   # TR=-0.1
    0.5603691749,   # TR=-0.01
    0.5737256299,   # TR=0
    0.5863493577,   # TR=+0.01
    0.6768917031,   # TR=+0.1
], dtype=float)

# CP-odd: sigma(TI) points, TR=0
TI_pts = np.array([-0.1, -0.01, 0.0, +0.01, +0.1], dtype=float)
sig_TI = np.array([
    0.5921394745,   # TI=-0.1
    0.5737580069,   # TI=-0.01
    0.5737256299,   # TI=0
    0.5739099319,   # TI=+0.01
    0.5919500347,   # TI=+0.1
], dtype=float)


# ============================================================
# Luminosity: 1 ab^-1 = 1e6 pb^-1
# ============================================================
L_pb = 1e6

# Asimov "data" = SM expectation
B = L_pb * sigma_SM


def q_asimov(mu, B):
    """Asimov -2 ln L ratio for Poisson, observed n=B, prediction mu."""
    if mu <= 0:
        return float("inf")
    return 2.0 * (mu - B - B * math.log(mu / B))


# ============================================================
# Fits
#   sigma(TR) = s0 + a TR + b TR^2
#   sigma(TI) = s0 + c TI^2
# ============================================================
X_TR = np.column_stack([np.ones_like(TR_pts), TR_pts, TR_pts**2])
coef_TR, *_ = np.linalg.lstsq(X_TR, sig_TR, rcond=None)
s0_TR, a_TR, b_TR = coef_TR


def sigma_of_TR(TR):
    return s0_TR + a_TR * TR + b_TR * TR**2


X_TI = np.column_stack([np.ones_like(TI_pts), TI_pts**2])
coef_TI, *_ = np.linalg.lstsq(X_TI, sig_TI, rcond=None)
s0_TI, c_TI = coef_TI


def sigma_of_TI(TI):
    return s0_TI + c_TI * TI**2


# ============================================================
# q(TR) and q(TI)
# ============================================================
def q_TR(TR):
    mu = L_pb * sigma_of_TR(TR)
    return q_asimov(mu, B)


def q_TI(TI):
    mu = L_pb * sigma_of_TI(TI)
    return q_asimov(mu, B)


# ============================================================
# Root finder (bisection): q(theta)=3.84 (95% CL, 1 dof)
# ============================================================
def bisect_root(func, a, b, q_target=3.84, tol=1e-12, maxit=200):
    fa = func(a) - q_target
    fb = func(b) - q_target
    if fa * fb > 0:
        raise ValueError(f"No sign change in [{a},{b}] (fa={fa}, fb={fb})")
    for _ in range(maxit):
        m = 0.5 * (a + b)
        fm = func(m) - q_target
        if abs(fm) < 1e-12 or 0.5 * (b - a) < tol:
            return m
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def find_limits():
    q_target = 3.84

    # TR: bracket roots around 0
    TR_scan = np.linspace(-0.05, 0.05, 2001)
    vals = np.array([q_TR(x) - q_target for x in TR_scan])

    neg_bracket, pos_bracket = None, None
    for i in range(len(TR_scan) - 1):
        if vals[i] * vals[i + 1] < 0:
            a, b = TR_scan[i], TR_scan[i + 1]
            if b < 0:
                neg_bracket = (a, b)
            elif a > 0:
                pos_bracket = (a, b)

    if neg_bracket is None or pos_bracket is None:
        raise RuntimeError("Could not bracket TR limits. Expand TR_scan range.")

    TR_low = bisect_root(q_TR, neg_bracket[0], neg_bracket[1], q_target=q_target)
    TR_high = bisect_root(q_TR, pos_bracket[0], pos_bracket[1], q_target=q_target)

    # TI symmetric
    TI_scan = np.linspace(0.0, 0.2, 4001)
    vals = np.array([q_TI(x) - q_target for x in TI_scan])

    pos_bracket = None
    for i in range(len(TI_scan) - 1):
        if vals[i] * vals[i + 1] < 0:
            pos_bracket = (TI_scan[i], TI_scan[i + 1])
            break

    if pos_bracket is None:
        raise RuntimeError("Could not bracket TI limit. Expand TI_scan range.")

    TI_lim = bisect_root(q_TI, pos_bracket[0], pos_bracket[1], q_target=q_target)

    return TR_low, TR_high, TI_lim


def print_summary(TR_low, TR_high, TI_lim):
    def print_yield(label, sig):
        N = L_pb * sig
        dN = math.sqrt(N)
        print(f"{label:25s}  sigma={sig:.10f} pb   N={N:10.1f}   sqrt(N)={dN:8.1f}")

    print("\n=== LHmuC yields at L = 1 ab^-1 (stat-only) ===")
    print_yield("SM", sigma_SM)

    print("\n=== Fitted parameterizations (pb) ===")
    print(f"sigma(TR) = {s0_TR:.10f} + ({a_TR:.10f})*TR + ({b_TR:.10f})*TR^2")
    print(f"sigma(TI) = {s0_TI:.10f} + ({c_TI:.10f})*TI^2")

    print("\n=== Expected 95% CL limits (rate-only, stat-only, L=1 ab^-1) ===")
    print(f"TR in [{TR_low:+.6f}, {TR_high:+.6f}]")
    print(f"TI in [{-TI_lim:+.6f}, {TI_lim:+.6f}]")

    N_SM = B
    dN_SM = math.sqrt(N_SM)
    dsig_stat = dN_SM / L_pb
    print("\n=== SM statistical uncertainty at 1 ab^-1 ===")
    print(f"N_SM = {N_SM:.1f},  sqrt(N_SM) = {dN_SM:.1f}")
    print(f"delta_sigma_stat = {dsig_stat:.6e} pb")


# ============================================================
# Plotting helpers
# ============================================================
def savefig(fig, path_prefix, save):
    if not save:
        return
    path_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path_prefix) + ".pdf")
    fig.savefig(str(path_prefix) + ".png", dpi=200)
    print(f"[saved] {path_prefix}.pdf")
    print(f"[saved] {path_prefix}.png")


def plot_sigma_TR(TR_low, TR_high, out_prefix, save):
    TR_grid = np.linspace(-0.12, 0.12, 1000)
    sig_grid = np.array([sigma_of_TR(x) for x in TR_grid])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(TR_grid, sig_grid, label="fit: σ(TR)")
    ax.scatter(TR_pts, sig_TR, label="MG points", zorder=3)
    ax.axvspan(TR_low, TR_high, alpha=0.2, label="95% CL (rate-only)")

    ax.set_xlabel("TR")
    ax.set_ylabel("σ [pb]")
    ax.set_title(r"LH$\mu$C @ 1 ab$^{-1}$: σ(TR) fit and points")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    savefig(fig, out_prefix / "LHmuC_sigma_vs_TR", save)


def plot_sigma_TI(TI_lim, out_prefix, save):
    TI_grid = np.linspace(-0.12, 0.12, 1000)
    sig_grid = np.array([sigma_of_TI(x) for x in TI_grid])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(TI_grid, sig_grid, label="fit: σ(TI)")
    ax.scatter(TI_pts, sig_TI, label="MG points", zorder=3)
    ax.axvspan(-TI_lim, TI_lim, alpha=0.2, label="95% CL (rate-only)")

    ax.set_xlabel("TI")
    ax.set_ylabel("σ [pb]")
    ax.set_title(r"LH$\mu$C @ 1 ab$^{-1}$: σ(TI) fit and points")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    savefig(fig, out_prefix / "LHmuC_sigma_vs_TI", save)


def plot_q_curves(TR_low, TR_high, TI_lim, out_prefix, save):
    q_target = 3.84

    TR_grid = np.linspace(-0.01, 0.01, 800)
    qtr = np.array([q_TR(x) for x in TR_grid])

    TI_grid = np.linspace(0.0, 0.10, 800)
    qti = np.array([q_TI(x) for x in TI_grid])

    fig1, ax1 = plt.subplots(figsize=(8, 5))
    ax1.plot(TR_grid, qtr, label="q(TR)")
    ax1.axhline(q_target, linestyle="--", label="95% CL threshold (3.84)")
    ax1.axvline(TR_low, linestyle=":", label=f"TR low = {TR_low:+.4f}")
    ax1.axvline(TR_high, linestyle=":", label=f"TR high = {TR_high:+.4f}")
    ax1.set_xlabel("TR")
    ax1.set_ylabel("q")
    ax1.set_title(r"LH$\mu$C @ 1 ab$^{-1}$: q(TR)")
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=9)
    fig1.tight_layout()
    savefig(fig1, out_prefix / "LHmuC_q_vs_TR", save)

    fig2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot(TI_grid, qti, label="q(TI)")
    ax2.axhline(q_target, linestyle="--", label="95% CL threshold (3.84)")
    ax2.axvline(TI_lim, linestyle=":", label=f"TI = {TI_lim:+.4f}")
    ax2.set_xlabel("TI")
    ax2.set_ylabel("q")
    ax2.set_title(r"LH$\mu$C @ 1 ab$^{-1}$: q(TI)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=9)
    fig2.tight_layout()
    savefig(fig2, out_prefix / "LHmuC_q_vs_TI", save)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--save", action="store_true", help="Save plots to --outdir as PDF+PNG.")
    ap.add_argument("--outdir", default="plots_limits_LHmuC", help="Output directory for plots.")
    args = ap.parse_args()

    TR_low, TR_high, TI_lim = find_limits()
    print_summary(TR_low, TR_high, TI_lim)

    out_prefix = Path(args.outdir)
    plot_sigma_TR(TR_low, TR_high, out_prefix, args.save)
    plot_sigma_TI(TI_lim, out_prefix, args.save)
    plot_q_curves(TR_low, TR_high, TI_lim, out_prefix, args.save)

    plt.show()


if __name__ == "__main__":
    main()
