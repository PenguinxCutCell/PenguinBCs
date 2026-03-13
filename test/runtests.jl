using Test
using StaticArrays
using PenguinBCs

@testset "BorderConditions validation" begin
    bc = BorderConditions(; left=Periodic(), right=Periodic())
    @test validate_borderconditions!(bc, 1) === bc

    bad = BorderConditions(; left=Periodic(), right=Dirichlet(0.0))
    @test_throws ArgumentError validate_borderconditions!(bad, 1)
end

@testset "Inflow/Outflow in BorderConditions" begin
    bc = BorderConditions(; left=Periodic(), right=Periodic())
    @test validate_borderconditions!(bc, 1) === bc

    bc2 = BorderConditions(; left=Inflow(1.0), right=Outflow())
    @test validate_borderconditions!(bc2, 1) === bc2

    bad = BorderConditions(; left=Periodic(), right=Outflow())
    @test_throws ArgumentError validate_borderconditions!(bad, 1)
end

@testset "eval_bc" begin
    x = SVector(1.0, 2.0)
    @test eval_bc(3.5, x, 0.1) == 3.5
    @test eval_bc((x, y) -> x + y, x, 0.1) == 3.0
    @test eval_bc((x, y, t) -> x + y + t, x, 0.1) == 3.1
end

@testset "Traction / outlet boundary types" begin
    bc = BorderConditions(
        ; left=Traction(SVector(1.0, 2.0)),
        right=PressureOutlet(0.0),
        bottom=DoNothing(),
        top=Traction((x, y, t) -> SVector(x + t, y - t)),
    )
    @test validate_borderconditions!(bc, 2) === bc

    x = SVector(0.2, 0.3)
    @test eval_bc((bc.borders[:left]).value, x, 0.5) == SVector(1.0, 2.0)
    @test eval_bc((bc.borders[:top]).value, x, 0.5) == SVector(0.7, -0.2)
    @test eval_bc((PressureOutlet()).value, x, 0.0) == 0.0
end

@testset "Symmetry boundary type" begin
    bc = BorderConditions(; left=Dirichlet(0.0), right=Symmetry())
    @test bc.borders[:right] isa Symmetry
    @test bc.borders[:right] isa AbstractBoundary

    # Periodic validation remains unchanged.
    bad = BorderConditions(; left=Periodic(), right=Symmetry())
    @test_throws ArgumentError validate_borderconditions!(bad, 1)
end

@testset "InterfaceConditions" begin
    ic = InterfaceConditions(
        scalar=ScalarJump(1.0, 2.0, 0.5),
        flux=FluxJump(3.0, 4.0, 0.25),
    )
    @test ic.scalar isa ScalarJump
    @test ic.flux isa FluxJump
end

@testset "GibbsThomson construction and eval" begin
    gt = GibbsThomson(0.25; kinetic=0.05)
    @test gt.capillary == 0.25
    @test gt.kinetic == 0.05

    gt_cb = GibbsThomson((x, y, t) -> 0.1 + x; kinetic=(x, y, t) -> 0.02 + t)
    x = SVector(0.3, 0.7)
    t = 0.4
    @test eval_bc(gt_cb.capillary, x, t) ≈ 0.4
    @test eval_bc(gt_cb.kinetic, x, t) ≈ 0.42
end

@testset "AlloyEquilibrium construction" begin
    ic = AlloyEquilibrium(0.3, 1.2, -0.8)
    @test ic isa AlloyEquilibrium
    @test ic isa AbstractInterfaceBC
    @test ic.k_partition == 0.3
    @test ic.T_m == 1.2
    @test ic.m_liquidus == -0.8
end
