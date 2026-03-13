# API

## Boundary Types

```@docs
AbstractBoundary
Dirichlet
Neumann
PenguinBCs.Robin
Periodic
Traction
PressureOutlet
DoNothing
Inflow
Outflow
BorderConditions
validate_borderconditions!
```

## Interface Types

```@docs
AbstractInterfaceBC
ScalarJump
FluxJump
RobinJump
GibbsThomson
AlloyEquilibrium
InterfaceConditions
```

## Evaluation

```@docs
eval_bc
```
