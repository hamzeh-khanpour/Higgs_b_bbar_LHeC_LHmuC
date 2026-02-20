# Higgs_b_bbar_LHeC_LHmuC
Higgs_b_bbar_LHeC_LHmuC

# Higgs_b_bbar_LHeC_LHmuC

**Sensitivity comparison: LHeC vs LHμC to SMEFT (dim-6) operators in the Higgs–bottom sector** using charged-current Higgs production in lepton–proton collisions:
\[
\ell^- p \to h\,j\,\nu_\ell,\qquad h\to b\bar b,\quad \ell=e,\mu
\]

This repository provides an end-to-end, reproducible pipeline:
**MadGraph5_aMC@NLO → (MadSpin) → Pythia8 → Delphes → analysis → expected limits**  
with a baseline **“no extra backgrounds”** assumption (SM Higgs is the only expectation).

---
---

## Motivation
We aim to **compare the sensitivity** of:
- **LHeC** (electron–proton) and
- **LHμC** (muon–proton)

to probe **Standard Model Effective Field Theory (SMEFT)** effects in the **Higgs–bottom** interactions, focusing on **dimension-six** deformations of the effective \(h b\bar b\) vertex.

---

## Physics setup

### Signal channel (charged-current Higgs production)
\[
\ell^- p \to h\, j\, \nu_\ell,\qquad h \to b\bar b
\]
- Parton-level: \(j\in\{u,d,s,c,b,g\}\)

### EFT parameterization
Following the reference study, we deform the effective Yukawa-like interaction after EWSB with two parameters:
- \(T_R\): CP-even
- \(T_I\): CP-odd

A minimal working form:
\[
y_b^{\rm eff} = y_b \left[1 + 3(T_R + iT_I)\right]
\]
so the SM is recovered at \(T_R=T_I=0\).

> **Note:** In this baseline study, the EFT modifies the **rate** dominantly; kinematic shapes are expected to be nearly SM-like (acceptance changes should be small). We still provide optional shape-based limits for completeness.

---

## Collider configurations
We consider:
- **LHeC:** \(E_e=60\) GeV, \(E_p=7\) TeV  \(\Rightarrow \sqrt{s}\approx 1.30\) TeV
- **LHμC:** \(E_\mu=500\) GeV, \(E_p=7\) TeV \(\Rightarrow \sqrt{s}\approx 3.74\) TeV

Integrated luminosities reported by default:
- \( \mathcal{L} = \{1,\;2,\;10\}\ {\rm ab}^{-1}\)

---

## Analysis assumptions
For the **LHeC vs LHμC sensitivity comparison**, we use a simplified baseline:

- **No extra backgrounds** (ignore \(Zjj\nu\), top, \(\gamma^\*/g^\*\), etc.)
- The **SM** process is the **only expected “background”**:
  - Asimov expectation: \(n = N_{\rm SM}\)
  - EFT prediction: \(\mu(T_R,T_I) = N(T_R,T_I)\)
- Limits are based on:
  - **Rate-only (single-bin Poisson)** likelihood (default)
  - Optional **shape-based (binned Poisson)** using 1D histograms (e.g. \(M_{bb}\))
- Optional normalization systematics on \(N_{\rm SM}\):
  - none (statistics-only)
  - 5%
  - 10%

---

## Repository structure
Planned/standard layout (create as you go):

