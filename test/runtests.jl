using PfaffianSystems
using DataStructures
using Bijections
using Test

@testset "PfaffianSystems.jl" begin
    # Write your tests here.
end

using PfaffianSystems: makeTestVarsAndIdeal

using Symbolics: @variables
# @testset "AsirWrapper.jl" begin
#     x, dx, var2do = genVars("x", 3)
#     @variables y[1:2]::Rational
#     os = OrderedSet((x[2], x[1]))

#     @test !isnothing(isAsirAvailable())
#     @test vec2str(x) == "x1,x2,x3"
#     @test vec2str(y) == "y1,y2"
#     @test vec2str(os) == "x2,x1"
#     @test isequal(asir_derivative(x[1]*sin(x[2]), x[2]), x[1]*cos(x[2]))
#     @test isequal(asir_derivative([x[1]*exp(x[2]), x[3]*x[2]^2], x[2]), [x[1]*exp(x[2]), 2*x[3]*x[2]])
#     @test isequal(asir_reduce((x[1]^2 - x[2]^2)/(x[1] + x[2])), x[1] - x[2])
#     @test isequal(asir_reduce([(x[1]^2 - x[2]^2)/(x[1] + x[2]), (sin(x[1])*x[2] + x[2]^2)/x[2]]), [x[1] - x[2], sin(x[1]) + x[2]])
#     a = (x[1]-1)^2*(x[2] - x[3])^3
#     @test isequal(asir_fctr(a), [-1=>1, x[1]-1=>2, x[3]-x[2]=>3])
# end

# @testset "DiffOps.jl" begin
#     x, dx, var2do = genVars("x", 3)
#     @test length(x) == length(dx) == length(var2do) == 3
#     @test map(s->haskey(var2do, s), x) |> all
#     @test map(s->haskey(inv(var2do), s), dx) |> all

#     y, dy, var2do = addVars("y", var2do)
#     @test length(var2do) == 4
#     @test haskey(var2do, y)
#     @test haskey(inv(var2do), dy)

#     @test apply_do(dx[1], x[1], var2do) == 1
#     @test apply_do(dx[2]^2 + 1, sin(x[2]), var2do) == 0
#     @test isequal(apply_do(x[2], x[1], var2do), x[1]*x[2])
#     @test isequal(apply_do(2.0, x[1], var2do), 2.0*x[1])
#     @test isequal(apply_do(3, x[1], var2do), 3*x[1])

#     @test apply_do(dx[1], x[1], var2do; use_asir=true) == 1
#     @test apply_do(dx[2]^2 + 1, sin(x[2]), var2do; use_asir=true) == 0
#     @test isequal(apply_do(x[2], x[1], var2do; use_asir=true), x[1]*x[2])
#     @test isequal(apply_do(2.0, x[1], var2do; use_asir=true), 2.0*x[1])
#     @test isequal(apply_do(3, x[1], var2do; use_asir=true), 3*x[1])

#     @test isequal(dmul(dx[2], 1, var2do), dx[2])
#     @test isequal(dmul(1, x[1], var2do), x[1])
#     @test isequal(dmul(dx[2], x[3]^2, var2do), dx[2]*x[3]^2)
#     @test isequal(dmul(dx[2], x[2]*dx[1], var2do), dx[1] + dx[1]*dx[2]*x[2])
#     @test isequal(dmul(dx[2], (x[2]^2 + 2*x[1])*dx[1], var2do), (x[2]^2 + 2*x[1])*dx[1]*dx[2] + 2*x[2]*dx[1])
# end

# @testset "DIdeals.jl" begin
#     x, dx, var2do, I = @test_nowarn makeTestVarsAndIdeal()
#     b = exp(-x[1]^2)*sin(x[2])*x[3]
#     y, dy, var2do = addVars("y", 2, var2do)
#     J = DIdeal([dx[1]^2 + 1, x[2]*dy[2] - 2], var2do)
#     @test isequal(stdmon!(I, OrderedSet(x)), [dx[2], 1])
#     @test isZeroDimensional(I)
#     @test !isZeroDimensional(J)
#     I3 = eliminationIdeal(I, x[1:2])
#     @test isequal(I3.gens, [x[3]*dx[3] - 1])
#     @test isequal(I3.v2d, Bijection(x[3], dx[3]))
#     @test_nowarn intersectionIdeal(I, J)
#     @test_nowarn integrationIdeal(I, x[1:1])
#     # @test_nowarn integrationIdeal(I, x[1:1])
#     # @test integrationIdeal(I, x[3:3]) |> isnothing
#     @test_nowarn restrictionIdeal(I, x[1:1])
#     @test isequal(apply_ideal(I, b), [0, 0, 0])
# end

# @testset "PfaffSys.jl" begin
#     x, dx, var2do, I = makeTestVarsAndIdeal()
#     pf = @test_nowarn PfaffianSystem(I)
#     @test isequal(get_vars(pf), x)
#     @test isequal(get_dvars(pf), dx)
#     funcAs, vars = @test_nowarn buildFuncA(pf)
#     x_bar = [0, 2, 1]
#     @test funcAs[1](x_bar) == [0 0; 0 0]
#     @test funcAs[2](x_bar) == [0 1; -1 0]
#     @test funcAs[3](x_bar) == [1 0; 0 1]
#     # @test_nowarn integrate(pf, [1 0; 0 1], [1, 1, 1], [3, 2, 1])
#     @test_nowarn integratePf(pf, [1 0; 0 1], Dict(x[1]=>1, x[2]=>1, x[3]=>1), Dict(x[1]=>3, x[2]=>2, x[3]=>1))
#     # @test_nowarn integrate(pf, [1 0; 0 1], cat([[cos(0.1*t), sin(0.1*t), 1+0.1*t] for t = 1:10]...; dims=2)) 
#     @test_nowarn integratePf(pf, [1 0; 0 1], 
#     Dict(
#         x[1]=>[cos(0.1*t) for t = 1:10], 
#         x[2]=>[sin(0.1*t) for t = 1:10],
#         x[3]=>[1+0.1*t for t = 1:10]
#         ) 
#     )
#     @test isequal(applyStdMons(pf, cos(x[2])), [cos(x[2]); -sin(x[2])])
#     @test isequal(denomLCM(pf), x[3])
# end

@testset "WeylAlgebra.jl" begin
    D, x, dx = weyl_algebra(QQ, ["x"])
    @test isequal(x*x, x^2)
    @test isequal(x*dx, x*dx)
    @test isequal(dx*x, x*dx+1)
    @test isequal(dx*dx, dx^2)
    @test isequal(x^2*x, x^3)
    @test isequal(x^2*dx, x^2*dx)
    @test isequal(dx^2*x, x*dx^2+2*dx)
    @test isequal(dx^2*dx,dx^3)
    @test isequal(x*x^2, x^3)
    @test isequal(x*dx^2, x*dx^2)
    @test isequal(dx*x^2, x^2*dx+2*x)
    @test isequal(dx*dx^2, dx^3)
    @test isequal(x^3*x, x^4)
    @test isequal(x^3*dx, x^3*dx)
    @test isequal(dx^3*x, x*dx^3+3*dx^2)
    @test isequal(dx^3*dx, dx^4)
    @test isequal(x^2*x^2, x^4)
    @test isequal(x^2*dx^2, x^2*dx^2)
    @test isequal(dx^2*x^2, x^2*dx^2+4*x*dx+2)
    @test isequal(dx^2*dx^2, dx^4)
    @test isequal(x*x^3, x^4)
    @test isequal(x*dx^3, x*dx^3)
    @test isequal(dx*x^3, x^3*dx+3*x^2)
    @test isequal(dx*dx^3, dx^4)
    @test isequal(x^4*x, x^5)
    @test isequal(x^4*dx, x^4*dx)
    @test isequal(dx^4*x, x*dx^4+4*dx^3)
    @test isequal(dx^4*dx, dx^5)
    @test isequal(x^3*x^2, x^5)
    @test isequal(x^3*dx^2, x^3*dx^2)
    @test isequal(dx^3*x^2, x^2*dx^3+6*x*dx^2+6*dx)
    @test isequal(dx^3*dx^2, dx^5)
    @test isequal(x^2*x^3, x^5)
    @test isequal(x^2*dx^3, x^2*dx^3)
    @test isequal(dx^2*x^3, x^3*dx^2+6*x^2*dx+6*x)
    @test isequal(dx^2*dx^3, dx^5)
    @test isequal(x*x^4, x^5)
    @test isequal(x*dx^4, x*dx^4)
    @test isequal(dx*x^4, x^4*dx+4*x^3)
    @test isequal(dx*dx^4, dx^5)

    @test isequa



end