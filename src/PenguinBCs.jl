module PenguinBCs

using StaticArrays

export AbstractBoundary, Dirichlet, Neumann, Robin, Periodic, Traction, PressureOutlet, DoNothing, Symmetry, Inflow, Outflow, BorderConditions, validate_borderconditions!
export AbstractInterfaceBC, ScalarJump, FluxJump, RobinJump, HarmonicAnisotropy, GibbsThomson, AlloyEquilibrium, InterfaceConditions
export eval_bc

"""
Abstract supertype for boundary conditions on the domain boundary.
"""
abstract type AbstractBoundary end

"""
Dirichlet boundary condition `u = value`.
"""
struct Dirichlet <: AbstractBoundary
    value::Union{Function,Float64}
end

"""
Neumann boundary condition `∂ₙu = value`.
"""
struct Neumann <: AbstractBoundary
    value::Union{Function,Float64}
end

"""
Robin boundary condition `α*u + β*∂ₙu = value`.
"""
struct Robin <: AbstractBoundary
    α::Union{Function,Float64}
    β::Union{Function,Float64}
    value::Union{Function,Float64}
end

"""
Periodic boundary marker.
"""
struct Periodic <: AbstractBoundary end

"""
Prescribed Cauchy traction boundary condition `σn = value`.
"""
struct Traction{T} <: AbstractBoundary
    value::T
end

"""
Pressure-outlet boundary condition with exterior pressure `p_out`.

For Stokes solvers this maps to `σn = -p_out n`.
"""
struct PressureOutlet{T} <: AbstractBoundary
    value::T
end
PressureOutlet() = PressureOutlet(0.0)

"""
Homogeneous traction boundary condition `σn = 0`.
"""
struct DoNothing <: AbstractBoundary end

"""
Symmetry / free-slip boundary condition for vector Stokes solvers.

Interpreted on a boundary with unit normal `n` as:
    u⋅n = 0
    (I - n⊗n) σ(u,p) n = 0
"""
struct Symmetry <: AbstractBoundary end

"""
Advection inflow boundary condition `u⋅n < 0`: impose transported scalar value.
"""
struct Inflow{T} <: AbstractBoundary
    value::T
end

"""
Advection outflow boundary marker `u⋅n ≥ 0`: no imposed scalar data.
"""
struct Outflow{T} <: AbstractBoundary end
Outflow(::Type{T}) where {T} = Outflow{T}()
Outflow() = Outflow{Float64}()

"""
Container for boundary conditions indexed by side symbols.
"""
struct BorderConditions
    borders::Dict{Symbol,AbstractBoundary}
end



"""
Construct [`BorderConditions`](@ref) from keyword pairs such as
`left=Dirichlet(0.0), right=Neumann(1.0)`.
"""
function BorderConditions(; kwargs...)
    d = Dict{Symbol,AbstractBoundary}()
    for (k, v) in kwargs
        v isa AbstractBoundary || throw(ArgumentError("border `$k` must be an AbstractBoundary"))
        d[Symbol(k)] = v
    end
    return BorderConditions(d)
end



function _side_pairs(N::Int)
    if N == 1
        return ((:left, :right),)
    elseif N == 2
        return ((:left, :right), (:bottom, :top))
    elseif N == 3
        return ((:left, :right), (:bottom, :top), (:backward, :forward))
    end
    throw(ArgumentError("unsupported dimension N=$N; expected 1, 2, or 3"))
end

function _validate_periodic_pairing!(borders, N)
    pairs = _side_pairs(N)
    for (lo, hi) in pairs
        blo = get(borders, lo, nothing)
        bhi = get(borders, hi, nothing)
        plo = blo isa Periodic
        phi = bhi isa Periodic
        if xor(plo, phi)
            throw(ArgumentError("periodic boundaries must be paired: both `$lo` and `$hi` must be Periodic"))
        end
    end
    return nothing
end

"""
Validate periodic side pairing for [`BorderConditions`](@ref) in dimension `N`.

Supported dimensions are `1`, `2`, and `3`. Opposite sides must both be periodic
or both non-periodic.
"""
function validate_borderconditions!(bc::BorderConditions, N)
    _validate_periodic_pairing!(bc.borders, N)
    return bc
end



"""
Abstract supertype for interface conditions across an internal interface.
"""
abstract type AbstractInterfaceBC end

"""
Scalar jump condition `α₁*u₁ - α₂*u₂ = value`.
"""
struct ScalarJump <: AbstractInterfaceBC
    α₁::Union{Function,Float64}
    α₂::Union{Function,Float64}
    value::Union{Function,Float64}
end

"""
Flux jump condition `β₁*∂ₙu₁ - β₂*∂ₙu₂ = value`.
"""
struct FluxJump <: AbstractInterfaceBC
    β₁::Union{Function,Float64}
    β₂::Union{Function,Float64}
    value::Union{Function,Float64}
end

"""
Robin-type jump condition `α*[u] + β*[∂ₙu] = value`.
"""
struct RobinJump <: AbstractInterfaceBC
    α::Union{Function,Float64}
    β::Union{Function,Float64}
    value::Union{Function,Float64}
end

"""
2D harmonic anisotropy descriptor for Gibbs-Thomson coefficients.

Stored parameters:
- `ϵ`: anisotropy amplitude (scalar or callback)
- `m`: symmetry fold (`m=4` by default)
- `θ0`: orientation offset angle (scalar or callback)
- `use_stiffness`: when `true`, Stefan may use stiffness correction
  `1 + ϵ(1-m^2)cos(m(θ-θ0))` for capillarity.
"""
struct HarmonicAnisotropy
    ϵ::Union{Function,Float64}
    m::Int
    θ0::Union{Function,Float64}
    use_stiffness::Bool
end

function HarmonicAnisotropy(
    ϵ;
    m::Int=4,
    θ0=0.0,
    use_stiffness::Bool=true,
)
    m > 0 || throw(ArgumentError("HarmonicAnisotropy requires m > 0"))
    return HarmonicAnisotropy(ϵ, m, θ0, use_stiffness)
end

"""
Gibbs-Thomson thermodynamic correction used by Stefan solvers.

Base isotropic coefficients:
- `capillary`: curvature coefficient `σ_κ`
- `kinetic`: kinetic coefficient `μ_Γ`

Optional 2D harmonic anisotropy descriptors:
- `capillary_anisotropy`
- `kinetic_anisotropy`

The complete interface law is applied by Stefan code as:
`TΓ = Tm - capillary_eff*κΓ - kinetic_eff*VΓ`.
"""
struct GibbsThomson <: AbstractInterfaceBC
    capillary::Union{Function,Float64}
    kinetic::Union{Function,Float64}
    capillary_anisotropy::Union{Nothing,HarmonicAnisotropy}
    kinetic_anisotropy::Union{Nothing,HarmonicAnisotropy}
end

GibbsThomson(capillary, kinetic) = GibbsThomson(capillary; kinetic=kinetic)

function GibbsThomson(
    capillary;
    kinetic=0.0,
    capillary_anisotropy::Union{Nothing,HarmonicAnisotropy}=nothing,
    kinetic_anisotropy::Union{Nothing,HarmonicAnisotropy}=nothing,
)
    return GibbsThomson(capillary, kinetic, capillary_anisotropy, kinetic_anisotropy)
end

"""
Alloy equilibrium interface descriptor for binary-alloy sharp-interface models.

Stored coefficients:
- `k_partition`: partition ratio `C_{s\\Gamma} = k\\,C_{l\\Gamma}`
- `T_m`: reference melting temperature
- `m_liquidus`: liquidus slope in `T_\\Gamma = T_m + m\\,C_{l\\Gamma}`
"""
struct AlloyEquilibrium{Tk,Tm} <: AbstractInterfaceBC
    k_partition::Tk
    T_m::Tm
    m_liquidus::Tm
end

"""
Container for scalar and flux interface conditions.
"""
struct InterfaceConditions
    scalar::Union{Nothing,AbstractInterfaceBC}
    flux::Union{Nothing,AbstractInterfaceBC}
end

"""
Construct [`InterfaceConditions`](@ref) from optional `scalar` and `flux` entries.
"""
InterfaceConditions(; scalar=nothing, flux=nothing) = InterfaceConditions(scalar, flux)

"""
Evaluate a boundary/interface value at spatial point `x` and time `t`.

`v` may be:
- `Float64` (returned directly),
- a callback `(x...)`,
- or a callback `(x..., t)`.
"""
eval_bc(v::Real, x::SVector, t) = Float64(v)
eval_bc(v::Ref, x::SVector, t) = eval_bc(v[], x, t)
eval_bc(v::AbstractVector, x::SVector, t) = v

function eval_bc(v::Function, x::SVector, t)
    y = nothing
    if applicable(v, x..., t)
        y = v(x..., t)
    elseif applicable(v, x...)
        y = v(x...)
    else
        throw(ArgumentError("boundary callback must accept (x...) or (x..., t)"))
    end
    return y isa Real ? Float64(y) : y
end

eval_bc(v, x::SVector, t) = eval_bc(Float64(v), x, t)

end
