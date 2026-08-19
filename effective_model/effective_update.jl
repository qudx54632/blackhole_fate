# Update one characteristic diamond in the effective model.
#
# The local equations for P and Y are defined in effective_functions.jl.

function effectiveDiamondStep(
    south,
    east,
    west,
    zMinusCenter,
    zPlusCenter,
    deltaMinus,
    deltaPlus;
    kwargs...,
)
    step = predictorCorrectorStep(
        south,
        east,
        west,
        deltaMinus,
        deltaPlus;
        kwargs...,
    )

    step.converged || error(
        "Fixed-point corrector failed at " *
        "(zMinus, zPlus) = ($zMinusCenter, $zPlusCenter).\n" *
        "iterations = $(step.iterations), residual = $(step.residual).",
    )

    ZNorth = step.north[1]
    xNorth = step.north[2]

    isfinite(ZNorth) || error(
        "Z is not finite at " *
        "(zMinus, zPlus) = ($zMinusCenter, $zPlusCenter).",
    )

    phiNorth = phiFromZ(ZNorth)
    phiNorth > 0.0 || error(
        "Phi underflowed at " *
        "(zMinus, zPlus) = ($zMinusCenter, $zPlusCenter).",
    )
    thetaNorth = thetaFromZX(ZNorth, xNorth)

    return (
        Z = ZNorth,
        phi = phiNorth,
        x = xNorth,
        theta = thetaNorth,
        north = step.north,
        iterations = step.iterations,
        residual = step.residual,
    )
end
