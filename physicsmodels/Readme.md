# Examples

Each script is self-contained, prints a validation table, writes CSV data and
saves figures to `examples/output/`.

```bash
# one-off setup (from the repository root): creates examples/Project.toml
julia --project=examples -e 'using Pkg; Pkg.develop(path="."); Pkg.add(["CairoMakie", "Printf"])'

# run any example
julia --project=examples examples/00_quickstart.jl
```

| script | what it shows | typical run time |
|:--|:--|:--|
| `quickstart.jl` | classify a mini-channel, friction, pressure drop, entrance lengths | seconds |
| `fully_developed_rectangular.jl` | exact vs polynomial vs finite-volume `fRe`, `Nu_H1`; 2nd-order convergence; velocity contours | ~10 s |
| `thermal_entrance.jl` | Graetz–Nusselt problem in rectangles, `Nu(x*)`, thermal entrance length | ~30 s |
| `developing_flow_2d.jl` | developing flow in a plane mini-channel, entrance length, `K∞` extrapolation | ~1–2 min |
| `heatsink_optimization.jl` | optimum channel width vs pressure budget | ~20 s |
| `tuckerman_pease.jl` | reproduction of Tuckerman & Pease (1981) Table I | seconds |
