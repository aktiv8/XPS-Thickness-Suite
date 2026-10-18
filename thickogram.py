"""
Thickogram: overlayer thickness from an overlayer peak and a substrate peak measured from DIFFERENT core levels.

Method: Peter J. Cumpson, Surf. Interface Anal. 29, 403-406 (2000). Python/Tk conversion of the GlowScript program
by Prof. Jong Wan Lee (Hallym University), then of the user's THickogram_v1_3.py (kept unmodified in legacy/).
See NOTICE for attribution and the licensing note.

Model. With x = t / (L_o cos(theta)) (theta from the surface NORMAL, 0 = normal emission), the substrate electrons
see the overlayer with an attenuation length that scales as E^0.75, so

    A = Io / (Is * Ro) = 2 * exp((B^0.75 - 0.5) x) * sinh(x / 2),     Ro = So / Ss,   B = KE_o / KE_s

and the thickness is t = x * L_o * cos(theta). x is found by solving f(x) = 0 with
f(x) = ln(A/2) - (B^0.75 - 0.5) x - ln sinh(x/2); f falls strictly from +inf to -inf, so there is always exactly
one root for A, B > 0. For B = 1 (equal kinetic energies) this is exactly x = ln(1 + A), i.e. the Strohmeier
equation used by the Calculator tab.

NOTE ON THE ANGLE: this module (and the Thickogram tab) use theta from the surface NORMAL, as Cumpson does. The
Calculator tab uses theta from the SURFACE. normal_from_surface() / surface_from_normal() convert.

Pure numpy / math, no Tk and no pyplot, so it is unit-testable and can draw on any matplotlib Axes.
"""
import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Nomogram constants, verbatim from THickogram_v1_3.py
# ---------------------------------------------------------------------------
B_INTERCEPT = 4.2          # the B scale is laid out at x_B = (4.2 - B) / 0.6
B_STEP = 0.6
AXP = np.arange(0.1, 6.1, 0.1)                                    # small dots on the green curve
AXD = [0.1, 0.2, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]                    # labelled dots on the green curve
AYP = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
       0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6,
       7, 8, 9, 10, 20, 30]                                       # ticks / lines of the A scale
AYD = [0.01, 0.1, 1.0, 10.0]                                      # labelled A values
BX = np.arange(0.4, 3.2, 0.2)                                     # the B scale lines

NOMOGRAM_B_RANGE = (0.4, 3.0)      # B values that exist on the printed scale
NOMOGRAM_X_RANGE = (0.1, 6.0)      # x values that exist on the green curve


# ---------------------------------------------------------------------------
# Maths
# ---------------------------------------------------------------------------

def sinh(x):
    """Same definition as the original script (works on scalars and arrays)."""
    return (np.exp(x) - np.exp(-x)) / 2


def func(x):
    """Original helper: x*((4.2 - 0.6x)^0.75 - 0.5), i.e. x_B*(B^0.75 - 0.5) at x_B = (4.2 - B)/0.6."""
    return x * ((B_INTERCEPT - B_STEP * x) ** 0.75 - 0.5)


def x_of_b(b):
    """Position of the B scale line for energy ratio B."""
    return (B_INTERCEPT - b) / B_STEP


def f_ab(a: float, b: float, x: float) -> float:
    """f(x) whose root is the solution (the original `fAB`)."""
    return math.log(a / 2.0) - (b ** 0.75 - 0.5) * x - math.log(math.sinh(x / 2.0))


def forward_a(x: float, b: float) -> float:
    """The intensity ratio A that a layer of reduced thickness x gives (inverse of solve_x)."""
    return 2.0 * math.exp((b ** 0.75 - 0.5) * x) * math.sinh(x / 2.0)


def solve_x(a: float, b: float) -> float:
    """
    Solve f(x) = 0 for x = t/(L cos(theta)) by bracketing + bisection (f is strictly decreasing, so the root is
    unique). Raises ValueError for non-finite or non-positive A or B instead of looping or returning NaN.
    """
    for name, v in (("A", a), ("B", b)):
        if not (isinstance(v, (int, float, np.floating, np.integer)) and math.isfinite(v)) or v <= 0:
            raise ValueError(f"{name} must be a positive finite number (got {v!r}).")
    lo = 1e-300
    if f_ab(a, b, lo) <= 0:
        raise ValueError("Intensity ratio A is too small to give a measurable thickness.")
    hi = 1.0
    while f_ab(a, b, hi) > 0:
        hi *= 2.0
        if hi > 2000.0:
            raise ValueError("Intensity ratio A is too large to solve (thickness beyond any physical range).")
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if f_ab(a, b, mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-14 * hi:          # relative, so tiny roots (A ~ 1e-10, x ~ 1e-10) are as precise as large ones
            break
    return 0.5 * (lo + hi)


def ratio_a(io: float, is_: float, so: float, ss: float) -> float:
    """A = Io / (Is * Ro) with Ro = So / Ss (all must be > 0)."""
    for name, v in (("Io", io), ("Is", is_), ("So", so), ("Ss", ss)):
        if not math.isfinite(v) or v <= 0:
            raise ValueError(f"{name} must be a positive number (got {v!r}).")
    return io / (is_ * (so / ss))


def energy_ratio_b(be_o: float, be_s: float, source_ev: float) -> float:
    """B = KE_overlayer / KE_substrate with KE = source energy - binding energy."""
    ke_o, ke_s = source_ev - be_o, source_ev - be_s
    if ke_o <= 0 or ke_s <= 0:
        raise ValueError(f"Kinetic energies must be positive (got {ke_o:g} eV and {ke_s:g} eV): "
                         "check the binding energies and the X-ray source energy.")
    return ke_o / ke_s


def normal_from_surface(theta_surface_deg: float) -> float:
    """Calculator convention (from the surface, 90 = normal emission) -> Thickogram convention (from the normal)."""
    return 90.0 - theta_surface_deg


def surface_from_normal(theta_normal_deg: float) -> float:
    return 90.0 - theta_normal_deg


@dataclass
class ThickogramResult:
    x: float                      # t / (L cos theta)
    t_nm: float
    a: float
    b: float
    ro: float
    ke_o: float
    ke_s: float
    lam_nm: float
    theta_normal_deg: float
    warnings: List[str] = field(default_factory=list)


def thickness(io: float, is_: float, so: float, ss: float, be_o: float, be_s: float, source_ev: float,
              lam_nm: float, theta_normal_deg: float) -> ThickogramResult:
    """The whole calculation. theta is from the surface NORMAL (0 = normal emission)."""
    if not math.isfinite(lam_nm) or lam_nm <= 0:
        raise ValueError("Lambda (attenuation length) must be a positive number of nm.")
    if not math.isfinite(theta_normal_deg) or not (0.0 <= theta_normal_deg < 90.0):
        raise ValueError("Theta (from the surface normal) must be from 0 up to (not including) 90 degrees.")
    a = ratio_a(io, is_, so, ss)
    b = energy_ratio_b(be_o, be_s, source_ev)
    x = solve_x(a, b)
    t = x * lam_nm * math.cos(math.radians(theta_normal_deg))
    warnings: List[str] = []
    if not (NOMOGRAM_B_RANGE[0] <= b <= NOMOGRAM_B_RANGE[1]):
        warnings.append(f"B = {b:.3f} is outside the printed nomogram scale ({NOMOGRAM_B_RANGE[0]}-"
                        f"{NOMOGRAM_B_RANGE[1]}); the thickness is still calculated, but the red line runs off the scale.")
    if not (NOMOGRAM_X_RANGE[0] <= x <= NOMOGRAM_X_RANGE[1]):
        warnings.append(f"x = {x:.3f} is outside the green curve's marked range ({NOMOGRAM_X_RANGE[0]}-"
                        f"{NOMOGRAM_X_RANGE[1]}).")
    return ThickogramResult(x=x, t_nm=t, a=a, b=b, ro=so / ss, ke_o=source_ev - be_o, ke_s=source_ev - be_s,
                            lam_nm=lam_nm, theta_normal_deg=theta_normal_deg, warnings=warnings)


# ---------------------------------------------------------------------------
# Drawing (the original plot_results, on a supplied Axes instead of a pyplot window)
# ---------------------------------------------------------------------------

def solution_line(a: float, b: float):
    """
    Endpoints of the red construction line, ((X0, Y0), (X1, Y1)).

    It runs from the A scale on the left axis, (0, ln(A/2)), to the B scale line at X = 2*x_B. Its crossing with
    the green curve y = ln sinh(X/4) is the answer, at X = 2x. (The original script ended the line at X = 2x, the
    thickness position, which draws the wrong slope - the 'redline plots correctly' to-do in its docstring.)
    """
    xb = x_of_b(b)
    y0 = math.log(a / 2.0)
    return (0.0, y0), (2.0 * xb, y0 - xb * (b ** 0.75 - 0.5))


def draw_nomogram(ax, result: Optional[ThickogramResult] = None, title: bool = True):
    """Draw the thickogram on `ax`. With a result, add the red solution line, the marker and the value."""
    ax.set_facecolor("white")

    x_vals = np.linspace(0.01, 6, 600)
    y_vals = np.log(sinh(x_vals / 2))
    ax.plot(2 * x_vals, y_vals, color="green", linewidth=1)

    for x in AXP:
        ax.plot(2 * x, np.log(sinh(x / 2)), "go", markersize=4)
    for x in AXD:
        ax.plot(2 * x, np.log(sinh(x / 2)), "go", markersize=7)

    for k in range(4, len(AYP) - 1):
        b_vals = np.arange(0.4, 3.1, 0.1)
        x_line = (B_INTERCEPT - b_vals) / B_STEP
        y_line = np.log(AYP[k] / 2) - func(x_line)
        ax.plot(2 * x_line, y_line, color="blue", linewidth=1)

    for b in BX:
        x = (B_INTERCEPT - b) / B_STEP
        y0 = np.log(AYP[4] / 2) - func(x)
        y1 = np.log(AYP[-2] / 2) - func(x)
        ax.plot([2 * x, 2 * x], [y0, y1], color="blue", linewidth=1)

    for y in AYP:
        ax.plot([0, 0.18], [np.log(y / 2)] * 2, color="black", linewidth=1)

    for y in AYD:
        ax.text(-0.4, np.log(y / 2), f"{y}", verticalalignment="center", color="black")

    for x in AXD:
        ax.text(2 * x + 0.15, np.log(sinh(x / 2)) - 0.25, f"{x}", color="green")

    for i in range(1, len(AYD)):
        xl = (B_INTERCEPT - 3) / B_STEP
        xr = (B_INTERCEPT - 0.4) / B_STEP
        yl = np.log(AYD[i] / 2) - func(xl)
        yr = np.log(AYD[i] / 2) - func(xr)
        ax.text(2 * xl - 0.4, yl, f"{AYD[i]:.1f}", color="blue")
        ax.text(2 * xr + 0.2, yr, f"{AYD[i]:.1f}", color="blue")

    for i in range(len(BX)):
        if i % 2 == 0:
            x = (B_INTERCEPT - BX[i]) / B_STEP
            y = np.log(AYP[4] / 2) - func(x)
            ax.text(2 * x + 0.05, y - 0.3, f"{BX[i]:.1f}", color="blue")

    if result is not None:
        x = result.x
        yk = math.log(math.sinh(x / 2))
        ax.plot([2 * x, 2 * x], [yk - 0.3, yk + 0.3], color="red", linewidth=2)
        (xa, ya), (xe, ye) = solution_line(result.a, result.b)
        ax.plot([xa, xe], [ya, ye], color="red", linewidth=2)
        ax.text(2 * x + 0.6, 2.4, f"x = {x:.4f}\nt = {result.t_nm:.2f} nm", color="red", fontweight="bold")

    if title:
        ax.set_title("Thickogram XPS Thickness Calculation")
    ax.set_xlabel("t/L*cos(theta) (2x)")
    ax.set_ylabel("log(Intensity Ratio)")
    ax.grid(True, linestyle="--", alpha=0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    return ax
