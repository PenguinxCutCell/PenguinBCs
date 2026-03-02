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

`value`, `α`, and `β` can be:
- A `Float64`
- A function of space: `(x...) -> ...`
- A function of space and time: `(x..., t) -> ...`

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

## Equations

A concise mapping to common PDE notation:

- Dirichlet:
  $u = g$
- Neumann:
  $\partial_n u = g$
- Robin:
  $ \alpha u + \beta\,\partial_n u = g$
- Scalar jump across interface \(\Gamma\):
  $\alpha_1 u_1 - \alpha_2 u_2 = g$
- Flux jump across interface \(\Gamma\):
  $\beta_1\,\partial_n u_1 - \beta_2\,\partial_n u_2 = g$

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

x = SVector(0.25, 0.75)
t = 0.1
val = eval_bc((x, y, t) -> x + y + t, x, t)

ic = InterfaceConditions(
    scalar = ScalarJump(1.0, 1.0, 0.0),
    flux   = FluxJump(1.0, 2.0, 0.0),
)
```

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
