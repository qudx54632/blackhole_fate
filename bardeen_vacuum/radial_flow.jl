# Radial-flow ODEs for a fixed zMinus.

using OrdinaryDiffEq

"""
    radialFlowODE!(du, u, , unusedParameters, unusedZPlus)

Evolve the radius and its sensitivity:

    dr/dzPlus = alpha * F(r) / 2
    dJ/dzPlus = alpha * Fp(r) * J / 2
"""
function radialFlowODE!(du, u, unusedParameters, unusedZPlus)
    r = u[1]
    J = u[2]

    du[1] = (alpha / 2.0) * F(r)
    du[2] = (alpha / 2.0) * Fp(r) * J

    return nothing
end

"""
    solveRadialFlow(r0)

Solve one radial flow with r(zPlusMin) = r0 and J(zPlusMin) = 1.

The returned solution keeps the ODE solver's dense interpolation.  Values
at plotting or PDE grid points should be obtained by calling the solution
at those coordinates rather than by indexing its internal save points.
"""
function solveRadialFlow(r0)
    initialValues = [Float64(r0), 1.0]
    zPlusSpan = (zPlusMin, zPlusMax)
    problem = ODEProblem(radialFlowODE!, initialValues, zPlusSpan)

    return solve(
        problem,
        Tsit5();
        dense = true,
        abstol = 1.0e-10,
        reltol = 1.0e-10,
    )
end
