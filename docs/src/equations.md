# Equations

This page maps package types to common PDE boundary/interface formulations.

## Boundary Conditions

Let `u` be the unknown scalar field, `n` the outward unit normal, and `g` a prescribed function.

Dirichlet (`Dirichlet(g)`):
```math
u = g \quad \text{on } \partial\Omega_D
```

Neumann (`Neumann(g)`):
```math
\partial_n u = g \quad \text{on } \partial\Omega_N
```

Robin (`Robin(\alpha, \beta, g)`):
```math
\alpha u + \beta\,\partial_n u = g \quad \text{on } \partial\Omega_R
```

Periodic (`Periodic()`): opposite sides are constrained by periodic identification. `validate_borderconditions!(bc, N)` enforces that periodic sides are paired.

Traction (`Traction(\tau)`):
```math
\sigma n = \tau \quad \text{on } \partial\Omega_T
```

Pressure outlet (`PressureOutlet(p_{\mathrm{out}})`): for Stokes solvers this maps to
```math
\sigma n = -p_{\mathrm{out}}\,n \quad \text{on } \partial\Omega_{\mathrm{out}}
```

Do-nothing (`DoNothing()`):
```math
\sigma n = 0 \quad \text{on } \partial\Omega_0
```

Advection inflow (`Inflow(g)`, where `u\cdot n < 0`): impose transported scalar value
```math
\phi = g \quad \text{on inflow boundary}
```

Advection outflow (`Outflow()`, where `u\cdot n \ge 0`): no scalar data is imposed at the boundary.

## Interface Conditions

Let `\Gamma` be an interior interface with traces `(\cdot)_1` and `(\cdot)_2` on each side.

Scalar jump (`ScalarJump(\alpha_1, \alpha_2, g)`):
```math
\alpha_1 u_1 - \alpha_2 u_2 = g \quad \text{on } \Gamma
```

Flux jump (`FluxJump(\beta_1, \beta_2, g)`):
```math
\beta_1\,\partial_n u_1 - \beta_2\,\partial_n u_2 = g \quad \text{on } \Gamma
```

Robin-type jump (`RobinJump(\alpha, \beta, g)`):
```math
\alpha [u] + \beta [\partial_n u] = g \quad \text{on } \Gamma
```
where `[q] := q_1 - q_2`.

Binary-alloy equilibrium descriptor (`AlloyEquilibrium(k, T_m, m)`):
```math
C_{s\Gamma} = k\,C_{l\Gamma}
```
```math
T_\Gamma = T_m + m\,C_{l\Gamma}
```

## Runtime Evaluation Semantics

For all coefficients/values, `eval_bc` supports:
- constant scalars
- callbacks with signature `(x...)`
- callbacks with signature `(x..., t)`

with `x` passed as `StaticArrays.SVector` and expanded into scalar coordinates.
