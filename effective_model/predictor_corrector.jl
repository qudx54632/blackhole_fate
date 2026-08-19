# Predictor-corrector evolution across one characteristic diamond.
#
#                 N = (zMinus[i+1], zPlus[j+1])
#                / \
#   W = (zMinus[i+1], zPlus[j])   E = (zMinus[i], zPlus[j+1])
#                \ /
#                 S = (zMinus[i], zPlus[j])
#
# Each point stores [Z, X], where Z = log(Phi).

# Reconstruct Phi, X, pMinus, and pPlus at the diamond center.
function diamondCenter(south, east, west, north, deltaMinus, deltaPlus)
    phiSouth = phiFromZ(south[1])
    phiEast = phiFromZ(east[1])
    phiWest = phiFromZ(west[1])
    phiNorth = phiFromZ(north[1])

    phi = (phiSouth + phiEast + phiWest + phiNorth) / 4.0
    x = (south[2] + east[2] + west[2] + north[2]) / 4.0
    pMinusValue =
        (phiWest + phiNorth - phiSouth - phiEast) / (2.0 * deltaMinus)
    pPlusValue =
        (phiEast + phiNorth - phiSouth - phiWest) / (2.0 * deltaPlus)

    return (; phi, x, pMinus = pMinusValue, pPlus = pPlusValue)
end


# Given the fields at S, E, and W, calculate the fields at N.
function predictorCorrectorStep(
    south,
    east,
    west,
    deltaMinus,
    deltaPlus;
    absoluteTolerance = 1.0e-10,
    relativeTolerance = 1.0e-8,
    maxIterations = 20,
)
    deltaMinus > 0.0 || error("deltaMinus must be positive.")
    deltaPlus > 0.0 || error("deltaPlus must be positive.")

    # First guess: neglect the mixed derivatives P and Y.
    predictedPhi = phiFromZ(east[1]) +
                   phiFromZ(west[1]) -
                   phiFromZ(south[1])
    north = if predictedPhi > 0.0
        [Z(predictedPhi), east[2] + west[2] - south[2]]
    else
        east .+ west .- south
    end
    residual = Inf

    for iteration in 1:maxIterations
        centerData = diamondCenter(
            south,
            east,
            west,
            north,
            deltaMinus,
            deltaPlus,
        )

        # Solve the effective equations for P and Y at the center.
        mixed = solvePY(
            centerData.phi,
            centerData.x,
            centerData.pMinus,
            centerData.pPlus,
        )

        phiNorth =
            phiFromZ(east[1]) +
            phiFromZ(west[1]) -
            phiFromZ(south[1]) +
            deltaMinus * deltaPlus * mixed.P

        if !isfinite(phiNorth) || phiNorth <= 0.0
            return (
                north = north,
                converged = false,
                iterations = iteration,
                residual = Inf,
            )
        end

        xNorth = east[2] + west[2] - south[2] +
                 deltaMinus * deltaPlus * mixed.Y
        newNorth = [Z(phiNorth), xNorth]

        all(isfinite, newNorth) ||
            error("The corrected north value is not finite.")

        residual = maximum(abs, newNorth .- north)
        tolerance = absoluteTolerance +
                    relativeTolerance * max(1.0, maximum(abs, newNorth))

        if residual <= tolerance
            return (
                north = newNorth,
                converged = true,
                iterations = iteration,
                residual = residual,
            )
        end

        north = newNorth
    end

    return (
        north = north,
        converged = false,
        iterations = maxIterations,
        residual = residual,
    )
end
