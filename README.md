# Sensorless PMSM Drive Simulation

Torque-controlled PMSM driving a quadratic fan load (11 kW @ 1500 rpm), built up in three
stages as Jupyter notebooks: plant + current control with real feedback, then a back-EMF
observer running alongside it, then the observer closing the control loop (fully sensorless).

**Start with `THEORY.md`** (or `THEORY.pdf`) for the theory, the observer structure, how the
gains were computed, and the full parameter reference — everything needed to understand what
the notebooks below implement and why.

## Notebooks (read/run in this order)

1. `01_pmsm_machine_and_controller.ipynb` — dq plant + dq PI current control, real θ/ω feedback.
2. `02_pmsm_with_parallel_observer.ipynb` — adds the αβ Structure-B back-EMF observer + PLL,
   running passively alongside (validates observer tuning before it's trusted).
3. `03_pmsm_sensorless_control.ipynb` — controller feedback switched to the observer's θ̂/ω̂
   (fully sensorless); shows the drive working, and where/why it breaks (a speed reversal
   triggers a genuine θ+π lock the drive cannot recover from on its own — see `THEORY.md` §7).

## Files

- `THEORY.md` / `THEORY.pdf` — theory, observer design, gain derivations, parameter reference.
- `pmsm_common.py` — shared machine/controller/observer parameters, gain formulas, coordinate
  transforms, saturation/anti-windup, and the convergence-aware torque-schedule simulation
  driver used by all three notebooks. This is the single source of truth for every constant
  and gain quoted in `THEORY.md` — read it alongside the notebook cells that call into it.
- `requirements.txt` — Python package versions used.

## Setup

```
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python -m ipykernel install --user --name=pmsm-sensorless --display-name="PMSM Sensorless (.venv)"
```

Then open a notebook, select the "PMSM Sensorless (.venv)" kernel, and run all cells. Or from
the command line:

```
jupyter nbconvert --to notebook --execute --inplace 01_pmsm_machine_and_controller.ipynb
```

Notebook `03` includes one intentionally slow cell (a couple of minutes) — see its markdown
for why: the `-20 N·m` schedule step drives the sensorless observer into a genuine chaotic
lock, which is the point of that notebook, not a bug.
