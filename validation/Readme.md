# Validation

* `python_crosscheck/` – independent NumPy/SciPy/Numba implementation of every numerical method in the
  package. It shares no code with the Julia sources; the Julia package was ported from it algorithm by
  algorithm. `python crosscheck.py --full` regenerates all numbers quoted in the README.
* `results/python_reference_output.txt` – saved output of that script.
* `data/literature_values.csv` – benchmark values with sources.

The Julia test-suite (`test/`) asserts the same cases with explicit tolerances.

## Known deviations (documented, not tuned away)
| item | deviation | reason |
|:--|:--|:--|
| T&P sample 3, ΔP | +24 % | bulk-mean viscosity; near-wall fluid hotter/thinner |
| T&P sample 3, R | −10.5 % | 1-D fin model, 4-wall `Nu_H1`, assumed substrate thickness |
| plane-channel K∞ | +2.5 % (0.691 vs 0.674) | two-point `1/Re` extrapolation on moderate grids; K converges slowly (inlet corner singularity) |
| entrance length | +3–5 % vs Chen | 1 % threshold on a discrete centre-line profile |
