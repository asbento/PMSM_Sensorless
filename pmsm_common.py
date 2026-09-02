"""
Shared constants and small helpers for the PMSM sensorless-control notebook series.

This module intentionally holds only:
  - machine / load / controller-gain PARAMETERS (and the formulas used to derive them),
  - generic coordinate-transform, saturation and plotting helpers,
  - a convergence-aware "run this torque schedule" simulation driver.

The actual physics (plant ODE), controller and observer right-hand-sides are written
directly in each notebook's cells -- that's the point of building this up in stages.
See PLAN.md / PMSM_Sensorless_Model_Spec.pdf for the full derivations.

Fixed state-vector layout used by every notebook (so the driver/plot helpers below
can index into it without per-notebook special-casing):

    0: isd        1: isq        2: omega_mech   3: theta_r      4: xi_d        5: xi_q
    6: ia_hat      7: ib_hat     8: ea_hat       9: eb_hat      10: z_a        11: z_b
    12: theta_pll  13: omega_pll

Notebook 1 only uses states 0-5. Notebooks 2 and 3 use all 14.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp

# --------------------------------------------------------------------------------------
# Machine parameters (Literature/RL_Parametros_PMSM.png)
# --------------------------------------------------------------------------------------
Rs = 18e-3          # stator resistance, ohm/phase
Ld = 180e-6          # d-axis synchronous inductance, H
Lq = 176e-6          # q-axis synchronous inductance, H
psi_m = 0.053        # permanent-magnet flux linkage, Wb (per PLAN.md convention: this is
                     # the "raw" per-phase value; the power-invariant factor sqrt(3/2) is
                     # applied explicitly wherever psi_m couples into the dq/alphabeta
                     # equations -- see psi_m_pi below)
p = 10               # pole pairs

# Power-invariant flux constant (see PLAN.md "Modeling convention"): this is the value
# that actually appears in the torque equation, the q-axis back-EMF feedforward, and the
# alpha-beta ground-truth back-EMF -- NOT the raw psi_m above.
psi_m_pi = np.sqrt(3 / 2) * psi_m   # ~= 0.06491 Wb

# --------------------------------------------------------------------------------------
# Mechanical / load parameters (not in the literature -- user-selected, see PLAN.md)
# --------------------------------------------------------------------------------------
J = 0.05             # rotor + load moment of inertia, kg*m^2
B = 0.0              # viscous friction, N*m/(rad/s) -- none; only the fan load damps

# Fan load: 11 kW at 1500 rpm mechanical defines the quadratic gain k in
#   T_load(omega_mech) = k * omega_mech * |omega_mech|
P_rated = 11e3                                  # W
n_rated_rpm = 1500.0                            # rpm
omega_mech_rated = n_rated_rpm * 2 * np.pi / 60  # rad/s  (~157.08)
T_rated = P_rated / omega_mech_rated             # N*m    (~70.02, matches the torque list)
k_fan = T_rated / omega_mech_rated ** 2          # N*m / (rad/s)^2  (~2.836e-3)


def load_torque(omega_mech: float) -> float:
    """Quadratic fan load torque, opposing rotation in either direction, continuous
    through zero (T_load(0)=0, dT_load/domega(0)=0)."""
    return k_fan * omega_mech * abs(omega_mech)


# --------------------------------------------------------------------------------------
# Current-loop PI gains (user-specified technical-optimum tuning against tau_u)
# --------------------------------------------------------------------------------------
tau_u = 10e-6        # s, inverter-response time constant used only for gain design


def current_pi_gains(Ld_: float, Lq_: float, Rs_: float, tau_u_: float):
    """Kp = L/(2*tau_u), Ki = R/(2*tau_u), separately per axis."""
    Kp_d = Ld_ / (2 * tau_u_)
    Ki_d = Rs_ / (2 * tau_u_)
    Kp_q = Lq_ / (2 * tau_u_)
    Ki_q = Rs_ / (2 * tau_u_)
    return Kp_d, Ki_d, Kp_q, Ki_q


Kp_d, Ki_d, Kp_q, Ki_q = current_pi_gains(Ld, Lq, Rs, tau_u)

# --------------------------------------------------------------------------------------
# Voltage saturation
# --------------------------------------------------------------------------------------
Vmax = 600.0         # V, circular limit on |[usd,usq]| (and on the observer's |[va,vb]|)


def saturate_voltage(ud: float, uq: float, vmax: float = Vmax):
    """Circular voltage-vector saturation. Returns (ud_sat, uq_sat, is_saturated)."""
    mag = np.hypot(ud, uq)
    if mag > vmax:
        scale = vmax / mag
        return ud * scale, uq * scale, True
    return ud, uq, False


def antiwindup_term(u_unsat: float, u_sat: float, kp: float, tt_gain: float = 1.0):
    """Back-calculation anti-windup correction to ADD to an integrator's derivative:
    dxi/dt = Ki*error + antiwindup_term(...). Zero when not saturated (u_sat==u_unsat).
    tt_gain scales the tracking speed relative to 1/Kp (tt_gain=1 -> tracking time
    constant Tt = Kp/Ki-ish, a standard, unremarkable choice)."""
    return tt_gain * kp * (u_sat - u_unsat)


# --------------------------------------------------------------------------------------
# Structure-B (PI-corrected Luenberger) alpha-beta back-EMF observer gains
# --------------------------------------------------------------------------------------
L_obs = (Ld + Lq) / 2   # isotropic inductance approximation used by the observer
# Triple-pole placement bandwidth (see PLAN.md derivation). The observer tracks e_alpha,e_beta
# as if a slowly-varying disturbance, but they are actually ROTATING at omega_r (up to ~1571
# rad/s electrical at the 1500 rpm rated point) -- so omega_n_obs must sit well above the
# highest electrical frequency in the schedule, or the estimated back-EMF vector picks up a
# speed-proportional steady-state phase bias. Empirically (see notebook 02 discussion) this
# residual comes from the observer's isotropic-L_obs approximation of the real (slightly
# salient, Ld != Lq) machine, and shrinks roughly as 1/omega_n_obs without ever fully
# vanishing (a sinusoidal-disturbance-rejection limit, not a bug): 20000 -> ~0.16 rad,
# 40000 -> ~0.08 rad, 80000 -> ~0.04 rad at rated speed. 80000 rad/s is used as a practical
# balance (still well-conditioned for RK45), acknowledged as a continuous-time idealization.
omega_n_obs = 80000.0


def observer_gains(R_: float, L_: float, omega_n: float):
    """Triple real pole at s=-omega_n on the observer's linearized error dynamics.
    Returns (Ki_P, Ki_I, Ke_P, Ke_I)."""
    Ki_P = 3 * omega_n - R_ / L_
    Ke_P = omega_n ** 2 * L_
    Ki_I = 2 * omega_n ** 2
    Ke_I = omega_n ** 3 * L_
    return Ki_P, Ki_I, Ke_P, Ke_I


Ki_P, Ki_I, Ke_P, Ke_I = observer_gains(Rs, L_obs, omega_n_obs)

# --------------------------------------------------------------------------------------
# PLL gains (position/speed extraction from the observer's estimated back-EMF angle)
# --------------------------------------------------------------------------------------
# Unlike the observer (whose bandwidth vs. electrical frequency directly sets steady-state
# accuracy), the PLL's steady-state accuracy for a constant-speed input is independent of its
# bandwidth (a type-2 loop tracks a phase ramp with zero steady-state error regardless of
# omega_n, given time to settle) -- confirmed empirically: sweeping omega_n_pll from 500 to
# 15000 rad/s gave an IDENTICAL 0.0426 rad steady-state bias every time (set entirely by the
# observer). Worse, a high omega_n_pll is actively harmful: fed a large, sudden phase step
# (e.g. the theta+pi flip, notebook 03), the PLL's sin(.) phase detector is nonlinear, and a
# high loop gain overreacts to a large-signal step in a cycle-slip-like way -- empirically a
# sharp threshold near the rated electrical frequency (~1571 rad/s): negligible overshoot
# below ~1000-1200 rad/s, thousands of rad/s of spurious overshoot above ~1600 rad/s (15000
# rad/s gave an ~18,900 rad/s overshoot spike on a pure pi step). 800 rad/s keeps settling
# time (~9 ms) far faster than the ~100 ms mechanical dynamics, with comfortable margin below
# the instability threshold and zero cost to steady-state accuracy.
omega_n_pll = 800.0
zeta_pll = 0.707


def pll_gains(omega_n: float, zeta: float):
    Kp_pll = 2 * zeta * omega_n
    Ki_pll = omega_n ** 2
    return Kp_pll, Ki_pll


Kp_pll, Ki_pll = pll_gains(omega_n_pll, zeta_pll)

# --------------------------------------------------------------------------------------
# Coordinate transforms: dq <-> alpha-beta are a plain rotation by theta (electrical).
# (No Clarke transform anywhere in this project -- we never model discrete abc phases.)
# --------------------------------------------------------------------------------------
def dq_to_alphabeta(d, q, theta):
    c, s = np.cos(theta), np.sin(theta)
    alpha = d * c - q * s
    beta = d * s + q * c
    return alpha, beta


def alphabeta_to_dq(alpha, beta, theta):
    c, s = np.cos(theta), np.sin(theta)
    d = alpha * c + beta * s
    q = -alpha * s + beta * c
    return d, q


def wrap_to_pi(theta):
    """Wrap angle(s) to (-pi, pi], for plotting/error metrics only (never needed inside
    an ODE right-hand side, where theta is kept unwrapped)."""
    return (theta + np.pi) % (2 * np.pi) - np.pi


# --------------------------------------------------------------------------------------
# Parameter bundle passed into every ODE right-hand side / controller / observer.
# Tref / isq_ref / isd_ref are mutated by the scheduler between chunks.
# --------------------------------------------------------------------------------------
def make_params(**overrides) -> SimpleNamespace:
    p_ = SimpleNamespace(
        Rs=Rs, Ld=Ld, Lq=Lq, psi_m=psi_m, psi_m_pi=psi_m_pi, p=p,
        J=J, B=B, k_fan=k_fan,
        tau_u=tau_u, Kp_d=Kp_d, Ki_d=Ki_d, Kp_q=Kp_q, Ki_q=Ki_q,
        Vmax=Vmax,
        L_obs=L_obs, Ki_P=Ki_P, Ki_I=Ki_I, Ke_P=Ke_P, Ke_I=Ke_I,
        Kp_pll=Kp_pll, Ki_pll=Ki_pll,
        Tref=0.0, isq_ref=0.0, isd_ref=0.0,
    )
    for k, v in overrides.items():
        setattr(p_, k, v)
    return p_


def set_torque_reference(params: SimpleNamespace, Tref: float):
    """Update params.Tref/isq_ref/isd_ref for a new torque command (id*=0 strategy)."""
    params.Tref = Tref
    params.isq_ref = Tref / (params.p * params.psi_m_pi)
    params.isd_ref = 0.0


# --------------------------------------------------------------------------------------
# Convergence-aware torque-schedule simulation driver
# --------------------------------------------------------------------------------------
@dataclass
class ScheduleResult:
    t: np.ndarray                  # (N,)
    X: np.ndarray                  # (n_states, N)
    level_times: list = field(default_factory=list)   # time at which each level starts
    level_values: list = field(default_factory=list)  # the Tref for each level
    converged_flags: list = field(default_factory=list)  # whether each level converged in time


def simulate_schedule(rhs, x0, params, torque_levels,
                       dt_chunk: float = 5e-3, max_time: float = 2.0,
                       domega_dt_tol: float = 0.5, iq_tol: float = 0.05,
                       n_confirm: int = 2, verbose: bool = True, debug: bool = False,
                       rtol: float = 1e-8, atol: float = 1e-10,
                       max_wall_time_per_level: float = 120.0) -> ScheduleResult:
    """Run the plant/controller(/observer) `rhs(t, x, params)` through each torque level
    in `torque_levels`, integrating in dt_chunk-sized pieces (adaptive RK45 within each
    chunk) and advancing to the next level once both the mechanical speed has stopped
    changing (|domega/dt| < domega_dt_tol) and iq has settled onto its reference
    (|iq-iq_ref| < iq_tol) for n_confirm consecutive chunks -- or after max_time of
    SIMULATED time, whichever comes first (a warning is printed if that cutoff is hit).

    `max_wall_time_per_level` is a separate, real-time (wall-clock) safety valve,
    independent of `max_time`: some regimes (e.g. a chaotic near-zero-speed lock, see
    notebook 03) make the adaptive integrator take drastically smaller steps than
    elsewhere, so a `max_time` that is cheap everywhere else can still be far more
    expensive there than the simulated-time budget suggests. Whichever cutoff (simulated
    or wall-clock) is hit first ends that level's integration and marks it non-converged.
    """
    t_chunks = [np.array([0.0])]
    x_chunks = [x0.reshape(-1, 1)]
    t0 = 0.0
    x = np.array(x0, dtype=float)
    level_times, level_values, converged_flags = [], [], []

    for level in torque_levels:
        set_torque_reference(params, level)
        level_times.append(t0)
        level_values.append(level)

        ok_streak = 0
        elapsed = 0.0
        prev_omega = x[2]
        converged = False
        n_chunks = 0
        import time as _time
        level_wall0 = _time.time()
        wall_cutoff = False
        while elapsed < max_time:
            if _time.time() - level_wall0 > max_wall_time_per_level:
                wall_cutoff = True
                break
            chunk_wall0 = _time.time()
            sol = solve_ivp(rhs, (t0, t0 + dt_chunk), x, args=(params,),
                             method="RK45", rtol=rtol, atol=atol)
            t_chunks.append(sol.t[1:])
            x_chunks.append(sol.y[:, 1:])
            x = sol.y[:, -1]
            t0 = sol.t[-1]
            elapsed += dt_chunk
            n_chunks += 1
            if debug:
                print(f"    [debug] level={level} chunk={n_chunks} "
                      f"chunk_wall={_time.time()-chunk_wall0:.3f}s nfev={sol.nfev} "
                      f"level_wall={_time.time()-level_wall0:.1f}s omega={x[2]:.3f}")

            domega_dt = (x[2] - prev_omega) / dt_chunk
            prev_omega = x[2]
            diq = abs(x[1] - params.isq_ref)

            if abs(domega_dt) < domega_dt_tol and diq < iq_tol:
                ok_streak += 1
            else:
                ok_streak = 0
            if ok_streak >= n_confirm:
                converged = True
                break

        converged_flags.append(converged)
        if verbose and not converged:
            reason = (f"wall-clock cutoff ({max_wall_time_per_level}s real time, "
                      f"reached {elapsed:.3f}s of simulated time)" if wall_cutoff else
                      f"{max_time}s of simulated time")
            print(f"  [warning] level Tref={level} N*m did not converge within {reason} "
                  f"(last |domega/dt|={domega_dt:.3f}, |diq|={diq:.4f})")

    t_all = np.concatenate(t_chunks)
    X_all = np.concatenate(x_chunks, axis=1)
    return ScheduleResult(t=t_all, X=X_all, level_times=level_times,
                           level_values=level_values, converged_flags=converged_flags)


# --------------------------------------------------------------------------------------
# Plotting helper
# --------------------------------------------------------------------------------------
def plot_traces(t, series: dict, level_times=None, title: str = "", figsize=(11, 8)):
    """series: dict of {subplot_label: [(trace_label, y_array, linestyle), ...]}.
    Draws one stacked subplot per key, with a legend, and vertical dashed lines at
    each torque-step time if level_times is given."""
    import matplotlib.pyplot as plt

    n = len(series)
    fig, axes = plt.subplots(n, 1, sharex=True, figsize=figsize)
    if n == 1:
        axes = [axes]
    for ax, (label, traces) in zip(axes, series.items()):
        for trace in traces:
            trace_label, y, *style = trace
            ls = style[0] if style else "-"
            ax.plot(t, y, ls, label=trace_label, linewidth=1.2)
        if level_times:
            for lt in level_times:
                ax.axvline(lt, color="gray", linestyle=":", linewidth=0.8, alpha=0.6)
        ax.set_ylabel(label)
        ax.legend(loc="upper right", fontsize=8, ncol=min(4, sum(1 for _ in traces)))
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel("time [s]")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig, axes
