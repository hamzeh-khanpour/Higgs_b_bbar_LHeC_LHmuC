#!/usr/bin/env python3
# Example:
#   python Higgs_b_bbar_LHeC_LHmuC_plots_LHmuC.py --save --outdir plots_LHmuC --tag LHmuC_Emu500 --pt-max 600

import math
import gzip
import argparse
from pathlib import Path
import re

import numpy as np
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


# -----------------------------
# LHmuC: user inputs (your LHE paths)
# -----------------------------
SAMPLES = [
    ("SM",
     "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/LHmuC_Emu500_CChiggsbb_SM/Events/run_01/LHmuC_Emu500_CChiggsbb_SM.lhe"),
    ("EFT CP-even (TR≠0, TI=0)",
     "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/LHmuC_Emu500_CChiggsbb_EFT_CP_even_TR/Events/run_01/LHmuC_Emu500_CChiggsbb_EFT_CP_even_TR.lhe"),
    ("EFT CP-odd (TR=0, TI≠0)",
     "/home/hamzeh-khanpour/MG5_aMC_v3_6_6/LHmuC_Emu500_CChiggsbb_EFT_CP_odd_TI/Events/run_01/LHmuC_Emu500_CChiggsbb_EFT_CP_odd_TI.lhe"),
]

# Titles updated for mu+ p and Emu=500
TOP_TITLE_PT  = r"LH$\mu$C: $\mu^+ p \to \bar{\nu}_\mu\, h\, j$"
TOP_TITLE_MBB = r"LH$\mu$C: $\mu^+ p \to \bar{\nu}_\mu\, h\, j$"
SUBTITLE      = r"$h\to b\bar b$ (MG/LHE level, decay in ME)"


# -----------------------------
# Helpers: open text (gz or not)
# -----------------------------
def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")


# -----------------------------
# Parse xsec from <init> block (pb)
# -----------------------------
def read_lhe_xsec_pb(filepath):
    with open_text(filepath) as f:
        in_init = False
        init_lines = []
        for line in f:
            if "<init>" in line:
                in_init = True
                continue
            if "</init>" in line and in_init:
                break
            if in_init:
                s = line.strip()
                if s and not s.startswith("#"):
                    init_lines.append(s)

    if len(init_lines) < 2:
        raise RuntimeError(f"Could not parse <init> block in {filepath}")

    first = init_lines[0].split()
    nprup = int(first[-1])
    proc_lines = init_lines[1:1 + nprup]
    if len(proc_lines) < nprup:
        raise RuntimeError(f"Init block in {filepath} does not contain {nprup} process lines")

    xsec_total = 0.0
    for pl in proc_lines:
        parts = pl.split()
        xsec_total += float(parts[0])  # XSECUP
    return xsec_total


# -----------------------------
# Robust parse (TR, TI) from BLOCK ALPPARS in the pre-event part of the LHE
# -----------------------------
RE_BLOCK_ALPPARS = re.compile(r"^\s*BLOCK\s+ALPPARS\b", re.IGNORECASE)
RE_BLOCK_ANY     = re.compile(r"^\s*BLOCK\b", re.IGNORECASE)
RE_DECAY_ANY     = re.compile(r"^\s*DECAY\b", re.IGNORECASE)

def read_tr_ti_from_lhe(filepath):
    # Read only the part before the first <event>
    pre_lines = []
    with open_text(filepath) as f:
        for line in f:
            if line.lstrip().startswith("<event"):
                break
            pre_lines.append(line.rstrip("\n"))

    tr, ti = None, None
    in_alppars = False

    for line in pre_lines:
        s = line.strip()

        if RE_BLOCK_ALPPARS.match(s):
            in_alppars = True
            continue

        if not in_alppars:
            continue

        # stop when ALPPARS block ends
        if RE_BLOCK_ANY.match(s) or RE_DECAY_ANY.match(s) or s.startswith("<"):
            break

        if not s or s.startswith("#"):
            continue

        s_nocom = s.split("#", 1)[0].strip()
        parts = s_nocom.split()
        if len(parts) < 2:
            continue

        try:
            idx = int(parts[0])
            val = float(parts[1])
        except Exception:
            continue

        if idx == 1:
            tr = val
        elif idx == 2:
            ti = val

        if tr is not None and ti is not None:
            break

    return tr, ti


def fmt_param(x):
    if x is None:
        return "0"
    if abs(x) < 1e-15:
        return "0"
    return f"{x:.3g}"


def fmt_tr_ti(tr, ti):
    if tr is None:
        tr = 0.0
    if ti is None:
        ti = 0.0
    return tr, ti


# -----------------------------
# Event iterator (streaming)
# Returns: (event_weight, particles)
# -----------------------------
def iter_lhe_events(filepath):
    with open_text(filepath) as f:
        for line in f:
            if line.lstrip().startswith("<event"):
                header = next(f).strip()
                while header == "":
                    header = next(f).strip()

                h = header.split()
                nup = int(h[0])
                xwgtup = float(h[2])  # event weight

                particles = []
                for _ in range(nup):
                    pline = next(f).strip()
                    while pline == "":
                        pline = next(f).strip()
                    p = pline.split()
                    pid = int(p[0])
                    status = int(p[1])
                    m1 = int(p[2])
                    m2 = int(p[3])
                    px = float(p[6])
                    py = float(p[7])
                    pz = float(p[8])
                    E  = float(p[9])
                    particles.append({
                        "id": pid, "status": status, "m1": m1, "m2": m2,
                        "px": px, "py": py, "pz": pz, "E": E
                    })

                # consume until </event>
                for endline in f:
                    if endline.strip().startswith("</event>"):
                        break

                yield xwgtup, particles


def pt(px, py):
    return math.hypot(px, py)


def invmass(p1, p2):
    E  = p1["E"]  + p2["E"]
    px = p1["px"] + p2["px"]
    py = p1["py"] + p2["py"]
    pz = p1["pz"] + p2["pz"]
    m2 = E*E - (px*px + py*py + pz*pz)
    return math.sqrt(m2) if m2 > 0 else 0.0


# -----------------------------
# Extract observables:
# - leading b pT from Higgs decay
# - m_bb from Higgs decay
# -----------------------------
def extract_observables(filepath):
    xsec_pb = read_lhe_xsec_pb(filepath)
    tr, ti = read_tr_ti_from_lhe(filepath)
    tr, ti = fmt_tr_ti(tr, ti)

    ptb1_list, mbb_list, w_raw_list = [], [], []
    sum_w_raw = 0.0

    for w_evt, parts in iter_lhe_events(filepath):
        sum_w_raw += w_evt

        higgs_indices = [i for i, p in enumerate(parts, start=1) if p["id"] == 25]

        b_from_h = []
        for i, p in enumerate(parts, start=1):
            if p["status"] != 1:
                continue
            if abs(p["id"]) != 5:
                continue
            if (p["m1"] in higgs_indices) or (p["m2"] in higgs_indices):
                b_from_h.append(p)

        # Fallback (shouldn't happen with decay chain)
        if len(b_from_h) < 2:
            b_all = [p for p in parts if p["status"] == 1 and abs(p["id"]) == 5]
            if len(b_all) >= 2:
                b_from_h = b_all[:2]

        if len(b_from_h) < 2:
            continue

        b1, b2 = b_from_h[0], b_from_h[1]
        ptb1 = max(pt(b1["px"], b1["py"]), pt(b2["px"], b2["py"]))
        mbb  = invmass(b1, b2)

        ptb1_list.append(ptb1)
        mbb_list.append(mbb)
        w_raw_list.append(w_evt)

    ptb1_arr = np.asarray(ptb1_list, float)
    mbb_arr  = np.asarray(mbb_list, float)
    w_raw    = np.asarray(w_raw_list, float)

    scale = xsec_pb / sum_w_raw if sum_w_raw != 0 else 0.0
    info = {
        "xsec_pb": float(xsec_pb),
        "tr": float(tr),
        "ti": float(ti),
        "scale": float(scale),
        "events_used": int(len(w_raw)),
    }
    return ptb1_arr, mbb_arr, w_raw, info


# -----------------------------
# Histogram to dsigma/dx with stat errors (pb/GeV)
# -----------------------------
def hist_dsig(values, w_raw, scale, edges):
    values = np.asarray(values, float)
    w = np.asarray(w_raw, float) * scale
    H, _  = np.histogram(values, bins=edges, weights=w)
    H2, _ = np.histogram(values, bins=edges, weights=w*w)
    widths = np.diff(edges)
    y  = np.divide(H, widths, out=np.zeros_like(H, dtype=float), where=widths > 0)
    dy = np.divide(np.sqrt(H2), widths, out=np.zeros_like(H2, dtype=float), where=widths > 0)
    return y, dy


def ratio_with_err(num, den, dnum, dden):
    R = np.divide(num, den, out=np.full_like(num, np.nan, dtype=float), where=(den > 0))
    dR = np.full_like(R, np.nan, dtype=float)
    mask = (num > 0) & (den > 0)
    dR[mask] = R[mask] * np.sqrt((dnum[mask]/num[mask])**2 + (dden[mask]/den[mask])**2)
    return R, dR


def step_with_plateau(edges, y):
    x = np.r_[edges[:-1], edges[-1]]
    yy = np.r_[y, y[-1] if len(y) else 0.0]
    return x, yy


def annotate_panel(ax, lines, fontsize=16):
    txt = "\n".join(lines)
    ax.text(
        0.04, 0.32, txt, transform=ax.transAxes,
        va="top", ha="left", fontsize=fontsize,
        bbox=dict(fc=(1, 1, 1, 0), ec=(1, 1, 1, 1), linewidth=2)
    )


def plot_cms_like(
    out_prefix: Path,
    title: str,
    x_label: str,
    y_label: str,
    edges: np.ndarray,
    curves,
    sm_label="SM",
    logy=False,
    y_floor=1e-14,
    save=False
):
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(9.8, 8.8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.20)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    for a in (ax, rax):
        a.tick_params(which="both", direction="in", top=True, right=True)
        a.grid(True, which="both", alpha=0.3)

    sm = next((c for c in curves if c["label"] == sm_label), None)
    if sm is None:
        raise RuntimeError(f"SM curve '{sm_label}' not found")

    # --- Top panel
    for c in curves:
        y_plot = np.where(c["y"] > 0, c["y"], y_floor) if logy else c["y"]
        x_step, y_step = step_with_plateau(edges, y_plot)
        ax.step(
            x_step, y_step, where="post",
            label=c["label"],
            lw=c.get("lw", 2.8),
            ls=c.get("ls", "-"),
        )

    ax.set_title(title, fontsize=20)
    ax.set_ylabel(y_label, fontsize=18)
    if logy:
        ax.set_yscale("log")

    ax.legend(ncol=3 if len(curves) >= 3 else 2, fontsize=12, loc="upper right")
    annotate_panel(ax, [c["annot"] for c in curves], fontsize=16)

    # --- Ratio panel
    den_y, den_dy = sm["y"], sm["dy"]
    for c in curves:
        if c["label"] == sm_label:
            continue
        R, dR = ratio_with_err(c["y"], den_y, c["dy"], den_dy)
        x_step, r_step = step_with_plateau(edges, np.nan_to_num(R, nan=0.0))
        rax.step(x_step, r_step, where="post", lw=2.0, ls=c.get("ls", "-"), label=f"{c['label']}/{sm_label}")

        centers = 0.5*(edges[:-1] + edges[1:])
        widths = np.diff(edges)
        good = np.isfinite(R) & np.isfinite(dR)
        rax.bar(
            centers[good],
            2.0*dR[good],
            bottom=(R[good] - dR[good]),
            width=widths[good],
            alpha=0.22,
            align="center",
            edgecolor="none",
        )

    rax.axhline(1.0, lw=1.0, color="black")
    rax.set_xlabel(x_label, fontsize=18)
    rax.set_ylabel("EFT/SM", fontsize=18)
    rax.set_ylim(0.0, 1.8)
    rax.legend(fontsize=10, loc="best")

    fig.subplots_adjust(left=0.11, right=0.98, top=0.92, bottom=0.08, hspace=0.10)

    if save:
        out_prefix.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(out_prefix) + ".pdf")
        fig.savefig(str(out_prefix) + ".png", dpi=200)
        print(f"[saved] {out_prefix}.pdf")
        print(f"[saved] {out_prefix}.png")

    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--save", action="store_true", help="Save plots as PDF+PNG.")
    ap.add_argument("--outdir", default="plots_LHmuC", help="Output directory.")
    ap.add_argument("--tag", default="LHmuC_Emu500", help="Filename prefix/tag.")
    ap.add_argument("--pt-bins", type=int, default=60, help="Number of pT bins.")
    ap.add_argument("--pt-max", type=float, default=-1.0, help="Max pT (GeV). If <0, auto from data.")
    ap.add_argument("--mbb-bins", type=int, default=70, help="Number of Mbb bins.")
    ap.add_argument("--mbb-min", type=float, default=90.0, help="Min Mbb (GeV).")
    ap.add_argument("--mbb-max", type=float, default=160.0, help="Max Mbb (GeV).")
    args = ap.parse_args()

    outdir = Path(args.outdir)

    # ---- Read samples
    samples = []
    print("\n=== Reading samples ===")
    for label, path in SAMPLES:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Missing file: {path}")

        ptb1, mbb, w_raw, info = extract_observables(str(path))
        TR_s = fmt_param(info["tr"])
        TI_s = fmt_param(info["ti"])

        print(f"{label:28s}  xsec={info['xsec_pb']:.6g} pb  TR={TR_s}, TI={TI_s}  events_used={info['events_used']}")

        samples.append({
            "label": label,
            "ptb1": ptb1,
            "mbb": mbb,
            "w_raw": w_raw,
            "xsec_pb": info["xsec_pb"],
            "scale": info["scale"],
            "tr": info["tr"],
            "ti": info["ti"],
        })

    # ---- Build common binning
    all_pt = np.concatenate([s["ptb1"] for s in samples])
    if args.pt_max > 0:
        pt_max = args.pt_max
    else:
        pt_max = max(50.0, np.percentile(all_pt, 99.5) * 1.3)
    pt_edges = np.linspace(0.0, pt_max, args.pt_bins + 1)
    mbb_edges = np.linspace(args.mbb_min, args.mbb_max, args.mbb_bins + 1)

    # ---- Build curves
    curves_pt = []
    curves_mbb = []
    for s in samples:
        y_pt, dy_pt = hist_dsig(s["ptb1"], s["w_raw"], s["scale"], pt_edges)
        y_m,  dy_m  = hist_dsig(s["mbb"],  s["w_raw"], s["scale"], mbb_edges)

        is_sm = (s["label"] == "SM")
        ls = "--" if is_sm else "-"
        lw = 2.8

        TR_s = fmt_param(s["tr"])
        TI_s = fmt_param(s["ti"])
        annot = (rf"{s['label']}:  $\sigma_{{gen}}={s['xsec_pb']:.4g}\ \mathrm{{pb}}$, "
                 rf"$TR={TR_s}$, $TI={TI_s}$")

        curves_pt.append({"label": s["label"], "y": y_pt, "dy": dy_pt, "ls": ls, "lw": lw, "annot": annot})
        curves_mbb.append({"label": s["label"], "y": y_m,  "dy": dy_m,  "ls": ls, "lw": lw, "annot": annot})

    # ---- Plot pT with log-y
    title_pt = TOP_TITLE_PT + "\n" + SUBTITLE
    out_prefix_pt = outdir / f"{args.tag}_dsig_dptb1_CMSstyle"
    plot_cms_like(
        out_prefix=out_prefix_pt,
        title=title_pt,
        x_label=r"$p_T(b_1)$  [GeV]",
        y_label=r"$d\sigma/dp_T$  [pb / GeV]",
        edges=pt_edges,
        curves=curves_pt,
        sm_label="SM",
        logy=True,
        y_floor=1e-14,
        save=args.save
    )

    # ---- Plot Mbb (linear y)
    title_mbb = TOP_TITLE_MBB + "\n" + SUBTITLE
    out_prefix_mbb = outdir / f"{args.tag}_dsig_dMbb_CMSstyle"
    plot_cms_like(
        out_prefix=out_prefix_mbb,
        title=title_mbb,
        x_label=r"$M_{b\bar b}$  [GeV]",
        y_label=r"$d\sigma/dM_{b\bar b}$  [pb / GeV]",
        edges=mbb_edges,
        curves=curves_mbb,
        sm_label="SM",
        logy=False,
        save=args.save
    )

    plt.show()


if __name__ == "__main__":
    main()
