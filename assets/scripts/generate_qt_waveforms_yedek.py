#!/usr/bin/env python3
"""
generate_qt_waveforms.py
========================
Generate physiologically realistic pulsatile flow-rate waveforms for a 1 mm
diameter artery, scaled from large-vessel measurements via Murray cube law.

Two waveforms are produced:
  Option A - Biphasic  (cerebral/carotid territory, antegrade-only flow)
             Source: Holdsworth et al. (1999), common carotid artery (CCA)
             Ref:    Physiological Measurement 20(2):219-240
             DOI:    https://doi.org/10.1088/0967-3334/20/2/301

  Option B - Triphasic (peripheral/femoral territory, includes flow reversal)
             Source: Mikheev et al. (2020), superficial femoral artery (SFA)
             Ref:    J. Physics: Conf. Series 1683:022090
             DOI:    https://doi.org/10.1088/1742-6596/1683/2/022090

Murray cube law (xi=3.0, conservative steady-laminar value):
  Q_mean(1mm) = Q_mean(ref) * (r_1mm / r_ref)^3

Womersley number at 1 mm artery (r=0.5mm):
  alpha = r * sqrt(omega/nu) = 0.0005 * sqrt(2*pi*70/60 / 3.3e-6) ~= 0.745
  Since alpha << 2: quasi-Poiseuille, parabolic profile valid at all instants:
    u(r,t) = 2 * Q(t)/A * (1 - (r/r0)^2)

Outputs:
  assets/img/qt_biphasic_waveform.png
  assets/img/qt_triphasic_waveform.png
  Printed OpenFOAM table (50 points each) to stdout.

Usage:
  python3 assets/scripts/generate_qt_waveforms.py

Requirements:
  numpy, scipy, matplotlib (optional -- skipped if not installed)
"""

import os
import math
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib

# Output directories
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMG_DIR   = os.path.join(REPO_ROOT, "assets", "img")
os.makedirs(IMG_DIR,  exist_ok=True)

# Physical parameters (1 mm artery)
R0    = 0.0005       # vessel radius [m]
A     = math.pi * R0**2  # cross-section [m2]
nu    = 3.3e-6       # kinematic viscosity [m2/s]
T     = 0.857        # cardiac period [s] (70 bpm)
omega = 2.0 * math.pi / T

alpha = R0 * math.sqrt(omega / nu)
print(f"Womersley alpha = {alpha:.3f}  (< 2 -> quasi-Poiseuille, parabolic profile valid)")
print()

# Murray cube law scaling
r_CCA = 3.15e-3;  Q_mean_CCA = 6.0e-6   # m, m3/s (Holdsworth 1999)
r_SFA = 2.5e-3;   Q_mean_SFA = 2.07e-6  # m, m3/s (Mikheev 2020)

Q_mean_bip = Q_mean_CCA * (R0 / r_CCA)**3
U_mean_bip = Q_mean_bip / A

Q_mean_tri = Q_mean_SFA * (R0 / r_SFA)**3
U_mean_tri = Q_mean_tri / A

print("Murray cube law scaling (xi=3.0):")
print(f"  CCA -> 1mm: Q_mean = {Q_mean_bip*1e9:.3f} nL/s = {Q_mean_bip*1e6:.4f} uL/s  "      f"U_mean = {U_mean_bip*100:.4f} cm/s")
print(f"  SFA -> 1mm: Q_mean = {Q_mean_tri*1e9:.3f} nL/s = {Q_mean_tri*1e6:.4f} uL/s  "      f"U_mean = {U_mean_tri*100:.4f} cm/s")
print()

# Digitized waveform control points
# Option A: Holdsworth CCA biphasic waveform
# Digitized from Holdsworth et al. (1999), Table 4 + Figure 11.
# Original T_orig=0.917s (66 bpm) scaled to T=0.857s (70 bpm) via time-compression.
# Q_norm = Q(t)/Q_mean >= 0 always (biphasic: no flow reversal).
# Systolic peak: t=0.142s, Q_norm=3.93; dicrotic notch: t=0.371s, Q_norm=0.067
bip_t = np.array([
    0.000, 0.050, 0.085, 0.110, 0.130, 0.142, 0.165, 0.200, 0.250,
    0.300, 0.330, 0.345, 0.371, 0.395, 0.420, 0.450, 0.500, 0.580,
    0.670, 0.760, 0.857
])
bip_Qn = np.array([
    0.40, 0.42, 0.60, 1.20, 2.50, 3.93, 3.20, 2.00, 1.20,
    0.70, 0.50, 0.35, 0.067, 0.25, 0.55, 0.62, 0.58, 0.52,
    0.46, 0.43, 0.40   # last = first for periodic BC
])

# Option B: Mikheev SFA triphasic waveform
# Digitized from Mikheev et al. (2020), Figure 1.
# T=0.857s (70 bpm). TRIPHASIC: includes flow reversal (Q_norm < 0).
# Systolic peak: t=0.133s, Q_norm=12 (AU=12);
# Max reversal: t=0.285s, Q_norm=-3.2
tri_t = np.array([
    0.000, 0.040, 0.075, 0.100, 0.120, 0.133, 0.155, 0.180, 0.210,
    0.235, 0.260, 0.285, 0.310, 0.340, 0.365, 0.390, 0.420, 0.450,
    0.480, 0.520, 0.580, 0.650, 0.730, 0.800, 0.857
])
tri_Qn = np.array([
    0.15,  0.30,  2.00,  6.00, 10.00, 12.00,  9.00,  5.00,  1.50,
   -1.00, -2.80, -3.20, -2.50, -1.20, -0.30,  0.40,  1.20,  1.60,
    1.50,  1.20,  0.80,  0.50,  0.30,  0.18,  0.15   # last = first
])


def make_waveform(t_ctrl, Qn_ctrl, T_period, N_dense=200, N_foam=50):
    """
    Fit periodic CubicSpline, renormalize so trapezoidal mean=1.0, return arrays.

    Parameters
    ----------
    t_ctrl   : control-point times (first==last for periodicity)
    Qn_ctrl  : Q_norm at control points (first==last value)
    T_period : cardiac period [s]
    N_dense  : points for dense evaluation (saved as .npy)
    N_foam   : points for OpenFOAM table

    Returns
    -------
    t_dense, Qn_dense, t_foam, Qn_foam, cs_norm
    """
    # Periodic spline
    cs_raw = CubicSpline(t_ctrl, Qn_ctrl, bc_type="periodic")

    # Dense evaluation
    t_d = np.linspace(0, T_period, N_dense, endpoint=False)
    Qn_raw = cs_raw(t_d)

    # Renormalize to mean = 1.0
    mean_raw = np.trapezoid(Qn_raw, t_d) / T_period
    scale = 1.0 / mean_raw if abs(mean_raw) > 1e-12 else 1.0
    cs_norm = CubicSpline(t_ctrl, Qn_ctrl * scale, bc_type="periodic")
    Qn_d = cs_norm(t_d)
    mean_check = np.trapezoid(Qn_d, t_d) / T_period
    print(f"  Renorm factor = {scale:.4f},  mean(Q_norm) after = {mean_check:.5f}")

    # OpenFOAM table
    t_f  = np.linspace(0, T_period, N_foam, endpoint=False)
    Qn_f = cs_norm(t_f)

    return t_d, Qn_d, t_f, Qn_f, cs_norm


print("=" * 60)
print("Option A -- Holdsworth CCA biphasic")
print("=" * 60)
bip_td, bip_Qnd, bip_tf, bip_Qnf, bip_cs = make_waveform(bip_t, bip_Qn, T)
bip_Qd = bip_Qnd * Q_mean_bip
bip_Ud = bip_Qd / A
bip_Qf = bip_Qnf * Q_mean_bip
bip_Uf = bip_Qf / A

print()
print("=" * 60)
print("Option B -- Mikheev SFA triphasic")
print("=" * 60)
tri_td, tri_Qnd, tri_tf, tri_Qnf, tri_cs = make_waveform(tri_t, tri_Qn, T)
tri_Qd = tri_Qnd * Q_mean_tri
tri_Ud = tri_Qd / A
tri_Qf = tri_Qnf * Q_mean_tri
tri_Uf = tri_Qf / A


def print_foam_table(name, tf, Qnf, Qf, Uf):
    print()
    print("=" * 60)
    print(f"OpenFOAM table -- {name}")
    print("  {:>8}  {:>10}  {:>12}  {:>14}".format("t [s]", "Q_norm", "Q [m3/s]", "U_mean [m/s]"))
    print("  " + "-" * 50)
    for i in range(len(tf)):
        print(f"  {tf[i]:8.4f}  {Qnf[i]:10.4f}  {Qf[i]:12.4e}  {Uf[i]:14.6f}")


print_foam_table("Biphasic (Holdsworth CCA -> 1mm)",  bip_tf, bip_Qnf, bip_Qf, bip_Uf)
print_foam_table("Triphasic (Mikheev SFA -> 1mm)",    tri_tf, tri_Qnf, tri_Qf, tri_Uf)

# Matplotlib plots (optional)
try:
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def _plot(td, Qnd, Ud, t_ctrl, Qn_ctrl_orig, Q_mean, cs_norm, title, fname):
        """Two-panel figure: Q_norm(t) and physical Q(t)/U(t)."""
        fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        fig.suptitle(title, fontsize=12, fontweight="bold")

        # Panel 1: Q_norm
        ax = axes[0]
        ax.plot(td * 1000, Qnd, "C0-", lw=2, label="Q_norm(t) [spline]")
        ax.axhline(0, color="k", lw=0.5, ls="--")
        ax.axhline(1, color="gray", lw=0.8, ls=":", label="Q_norm = 1 (mean)")
        ax.set_ylabel("Q_norm = Q(t)/Q_mean", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.4)

        # Panel 2: velocity
        ax = axes[1]
        ax.plot(td * 1000, Ud * 100, "C1-", lw=2, label="U_mean(t) [cm/s]")
        ax.axhline(0, color="k", lw=0.5, ls="--")
        ax.set_xlabel("Time [ms]", fontsize=10)
        ax.set_ylabel("Mean velocity [cm/s]", fontsize=10)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.4)

        fig.tight_layout()
        out = os.path.join(IMG_DIR, fname)
        fig.savefig(out, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {out}")

    _plot(bip_td, bip_Qnd, bip_Ud, bip_t, bip_Qn, Q_mean_bip, bip_cs,
          "Biphasic waveform (Holdsworth CCA 1999 -> 1 mm artery)",
          "qt_biphasic_waveform.png")
    _plot(tri_td, tri_Qnd, tri_Ud, tri_t, tri_Qn, Q_mean_tri, tri_cs,
          "Triphasic waveform (Mikheev SFA 2020 -> 1 mm artery)",
          "qt_triphasic_waveform.png")

except ImportError:
    print()
    print("WARNING: matplotlib not installed -- skipping PNG output.")
    print("  Install with:  pip install matplotlib")
