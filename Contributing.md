# Contributing

Contributions that improve correctness, validation coverage or documentation are very welcome.

## Ground rules for a research code
1. **Every new model needs a validation case** — an exact solution, a manufactured solution, or published
   benchmark data with citation. Add it to `test/` and, if it produces a headline number, to
   `validation/` and the README table.
2. **State tolerances honestly.** If a model deviates from a benchmark, assert the deviation you observe
   and explain it in the docstring or docs — do not tune coefficients to hide it.
3. **No unverified coefficients.** Correlation constants must be cross-checked (numerically or against the
   primary source) before merging; cite the equation number.
4. Keep the core dependency-free (standard library only). Plotting belongs in `examples/`.

## Workflow
```bash
git clone https://github.com/start-again-06/MiniChannelFlow.jl && cd MiniChannelFlow.jl
julia --project -e 'using Pkg; Pkg.test()'
```
Format with [JuliaFormatter](https://github.com/domluna/JuliaFormatter.jl) (`.JuliaFormatter.toml` is provided):
`julia -e 'using JuliaFormatter; format(".")'`.

Open an issue before large changes. Pull requests should pass CI and update `CHANGELOG.md`.
