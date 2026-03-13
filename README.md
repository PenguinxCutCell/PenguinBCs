# PenguinBCs.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinBCs.jl/dev)
![CI](https://github.com/PenguinxCutCell/PenguinBCs.jl/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/PenguinxCutCell/PenguinBCs.jl/graph/badge.svg?token=Q50hXtjAKk)](https://codecov.io/gh/PenguinxCutCell/PenguinBCs.jl)


`PenguinBCs.jl` provides lightweight boundary-condition and interface-condition types for PDE solvers in the PenguinxCutCell ecosystem.

It focuses on two needs:
- A shared, explicit data model for boundary and interface constraints.
- A small runtime helper (`eval_bc`) that evaluates constant or callback-based values consistently.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/PenguinxCutCell/PenguinBCs.jl")
```

## Boundary Conditions

Supported boundary types:
- `Dirichlet(value)`
- `Neumann(value)`
- `Robin(α, β, value)`
- `Periodic()`
- `Traction(value)`
- `PressureOutlet(p_out)` or `PressureOutlet()`
- `DoNothing()`
- `Inflow(value)`
- `Outflow()`

`value`, `α`, and `β` can be:
- A `Float64`
- A function of space: `(x...) -> ...`
- A function of space and time: `(x..., t) -> ...`

Boundary-type intent:
- `Dirichlet`: prescribed state value.
- `Neumann`: prescribed normal derivative / scalar flux.
- `Robin`: mixed state/flux condition.
- `Periodic`: periodic side marker (must be paired with opposite side).
- `Traction`: prescribed Cauchy traction vector/tensor contribution.
- `PressureOutlet`: exterior pressure outlet (for Stokes, `σn = -p_out n`).
- `DoNothing`: homogeneous traction (`σn = 0`).
- `Inflow`: advection inflow marker with prescribed transported scalar.
- `Outflow`: advection outflow marker with no imposed scalar data.

`BorderConditions` stores boundary assignments by side symbol:
- 1D: `:left`, `:right`
- 2D: `:left`, `:right`, `:bottom`, `:top`
- 3D: `:left`, `:right`, `:bottom`, `:top`, `:backward`, `:forward`

Periodic boundaries must come in opposite-side pairs and can be validated with `validate_borderconditions!(bc, N)`.

## Interface Conditions

Supported interface jump conditions:
- `ScalarJump(α₁, α₂, value)`
- `FluxJump(β₁, β₂, value)`
- `RobinJump(α, β, value)`

`InterfaceConditions(; scalar=nothing, flux=nothing)` groups scalar and flux interface conditions.

## Esoteric BC's/IC's

- `GibbsThomson(capillary; kinetic=0.0)`
`GibbsThomson(capillary; kinetic=0.0)` stores isotropic Gibbs-Thomson correction terms:
- `capillary`: curvature coefficient
- `kinetic`: kinetic velocity coefficient
- Imposes a modified temperature on the interface from `Robin(1.0, 0.0, Tm)` with `TΓ = Tm - capillary*κΓ - kinetic*VΓ`

- `AlloyEquilibrium(k_partition, T_m, m_liquidus)`
`AlloyEquilibrium(k_partition, T_m, m_liquidus)` stores binary alloy equilibrium interface conditions:
- `k_partition`: partition coefficient \(C_{s\Gamma} = k\,C_{l\Gamma}\)
- `T_m`: melting temperature
- `m_liquidus`: liquidus slope \(T_\Gamma = T_m + m\,C_{l\Gamma}\)

## Equations

A concise mapping to common PDE notation:

- Dirichlet:
  $u = g$
- Neumann:
  $\partial_n u = g$
- Robin:
  $ \alpha u + \beta\,\partial_n u = g$
- Traction:
  $ \sigma n = \tau$
- Pressure outlet:
  $ \sigma n = -p_{\mathrm{out}} n$
- Do-nothing:
  $ \sigma n = 0$
- Advection inflow:
  $ \phi = g \text{ on } u\cdot n < 0$
- Advection outflow:
  no imposed scalar data on $u\cdot n \ge 0$
- Scalar jump across interface \(\Gamma\):
  $\alpha_1 u_1 - \alpha_2 u_2 = g$
- Flux jump across interface \(\Gamma\):
  $\beta_1\,\partial_n u_1 - \beta_2\,\partial_n u_2 = g$
- Binary-alloy equilibrium interface:
  $C_{s\Gamma} = k\,C_{l\Gamma}$
  $T_\Gamma = T_m + m\,C_{l\Gamma}$

## Quick Example

```julia
using StaticArrays
using PenguinBCs

bc = BorderConditions(
    left   = Dirichlet((x, y, t) -> sin(t)),
    right  = Neumann(0.0),
    bottom = Periodic(),
    top    = Periodic(),
)
validate_borderconditions!(bc, 2)

# Additional boundary families
bc_stokes = BorderConditions(
    left   = Traction(SVector(1.0, 2.0)),
    right  = PressureOutlet(0.0),
    bottom = DoNothing(),
    top    = Traction((x, y, t) -> SVector(x + t, y - t)),
)
validate_borderconditions!(bc_stokes, 2)

bc_adv = BorderConditions(
    left  = Inflow((x, t) -> 1 + x + t),
    right = Outflow(),
)
validate_borderconditions!(bc_adv, 1)

x = SVector(0.25, 0.75)
t = 0.1
val = eval_bc((x, y, t) -> x + y + t, x, t)
traction_val = eval_bc((bc_stokes.borders[:top]).value, x, t)

ic = InterfaceConditions(
    scalar = ScalarJump(1.0, 1.0, 0.0),
    flux   = FluxJump(1.0, 2.0, 0.0),
)

alloy_ic = AlloyEquilibrium(0.2, 0.0, -0.5)
```

## Runtime Evaluation (`eval_bc`)

`eval_bc(v, x, t)` supports:
- constants: `v::Real`
- spatial callbacks: `(x...) -> ...`
- space-time callbacks: `(x..., t) -> ...`
- vector-like values (returned as-is)

`x` is passed as `StaticArrays.SVector` and expanded to scalar coordinates when calling callback-based values.

## Documentation

Build docs locally:

```julia
julia --project=docs docs/make.jl
```

Generated pages include equations, references, and full API documentation.

## Development

Run tests:

```julia
julia --project -e 'using Pkg; Pkg.test()'
```
