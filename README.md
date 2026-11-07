# HALO — Handover-Anchored Location Optimisation for GSM

Simulator, experiment scripts and results for the manuscript

> K. Popat, D. Meva and H. Molia, "HALO: Handover-Anchored Location Optimisation for Joint Handover and Location Management in GSM Networks" (under submission).

HALO is a network-side scheme that keeps a per-user mobility record built from handover and location-update events and reuses it for
**ILR** (implicit location re-registration via TMSI reallocation), **HAP** (handover-anchored two-step paging) and
**DAH** (dwell-adaptive handover margin with a ping-pong return guard).

## Repository layout

| Path | Contents |
|---|---|
| `sim.py` | Time-stepped system-level GSM simulator: 144 omni cells on a hexagonal torus, COST-231 Hata (900 MHz), Gudmundson-correlated shadowing (σ = 8 dB, inter-site correlation 0.5), reuse-4 C/I, four mobility classes, Poisson MO/MT calls. All schemes run in lockstep on common random numbers. |
| `experiments.py` | All experiments in the paper (E1 main comparison, E2 LA-size sweep, E3 traffic intensity, E4 mobility mix, E5 paging delay penalty). Evaluation seeds 1–10. Resumable. |
| `analyze.py` | Means, 95 % Student-t intervals, paired Wilcoxon signed-rank tests, figures. |
| `make_numbers.py` | Derives every number quoted in the manuscript from `results/summary.json`. |
| `results/raw.jsonl` | Raw per-replication output (one JSON line per experiment × seed). |
| `results/summary.json`, `results/tests_E1.json` | Aggregated statistics and paired tests. |
| `figs/` | Figures 1–5 of the manuscript (300 dpi). |
| `tuning/` | Exploratory scripts used on the separate tuning seed 101 only (not used for any reported result). |

## Schemes

`GSM` (standard), `VAHO` (velocity-aware handover with oracle speed), `TSP` (two-step last-cell paging), `VAHO+TSP`, `HALO`, and ablations `HALO-noILR`, `HALO-noHAP`, `HALO-noDAH`.

## Reproducing the results

```bash
python -m pip install -r requirements.txt
python sim.py 1            # smoke test: 200 users, 25 simulated minutes (~10 s)
python experiments.py      # full set: 100 runs, ~35 min on 2 CPU cores (resumes if interrupted)
python analyze.py          # writes results/summary.json, results/tests_E1.json and figs/
python make_numbers.py     # writes numbers.json used by the manuscript
```

To regenerate from scratch, delete `results/raw.jsonl` first; otherwise `experiments.py` skips runs already present.
Results are deterministic for a given seed, NumPy version and platform; small floating-point differences across NumPy versions can shift individual counts slightly without changing the conclusions.

Tested with Python 3.13, NumPy 2.5, SciPy 1.18, Matplotlib 3.11 (Python ≥ 3.10 should work).

## Citation

See `CITATION.cff`. Please cite the paper if you use this code.

## License

MIT — see `LICENSE`.
