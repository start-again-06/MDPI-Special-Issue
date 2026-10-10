# Changelog

All notable changes are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/) and the project uses
[Semantic Versioning](https://semver.org/).

## [0.1.0] – 2026-10-08

### Added
- Fluid properties (`water(T)`), rectangular-channel geometry and Kandlikar–Grande size classification.
- Shah & London polynomial correlations (`fRe`, `Nu_T`, `Nu_H1`), Muzychka–Yovanovich apparent friction factor, Chen entrance length.
- Exact double-Fourier-series solution for fully developed rectangular-duct flow.
- Second-order finite-volume cross-section solver: `fRe`, `Nu_H1`.
- Thermal-entrance (Graetz–Nusselt) marching solver for rectangular ducts, constant wall temperature.
- 2-D incompressible Navier–Stokes solver (MAC projection method) for developing plane-channel flow.
- Parallel-channel heat-sink model, fixed-ΔP flow solver and channel-width optimiser.
- Test suite, examples, Documenter.jl site, independent Python cross-check implementation.
