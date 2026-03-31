# PenguinBCs.jl

`PenguinBCs.jl` defines boundary and interface condition types shared across Penguin operators and physics packages.

## Scope

The package provides:
- Boundary-condition types: `Dirichlet`, `Neumann`, `Robin`, `Periodic`, `Traction`, `PressureOutlet`, `DoNothing`, `Inflow`, `Outflow`
- Interface-condition types: `ScalarJump`, `FluxJump`, `RobinJump`, `AlloyEquilibrium`
- Gibbs-Thomson support: `GibbsThomson` with optional `HarmonicAnisotropy`
- Condition containers: `BorderConditions`, `InterfaceConditions`
- Evaluation helper: `eval_bc`

## Quick Start

```julia
using StaticArrays
using PenguinBCs

bc = BorderConditions(; left=Dirichlet(1.0), right=Neumann(0.0))
validate_borderconditions!(bc, 1)

x = SVector(0.5)
t = 0.0
val = eval_bc((x, t) -> x + t, x, t)
```

## Next

- Mathematical forms: [Equations](equations.md)
- Background material: [References](references.md)
- Full public API: [API](api.md)
