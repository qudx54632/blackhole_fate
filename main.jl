# Main entry point
using Serialization

include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "bardeen_vacuum", "functions.jl"))
include(joinpath(@__DIR__, "bardeen_vacuum", "horizons.jl"))
include(joinpath(@__DIR__, "grid.jl"))
include(joinpath(@__DIR__, "bardeen_vacuum", "radial_flow.jl"))
include(joinpath(@__DIR__, "bardeen_vacuum", "boundary_values.jl"))
include(joinpath(@__DIR__, "effective_model", "effective_functions.jl"))
include(joinpath(@__DIR__, "effective_model", "predictor_corrector.jl"))
include(joinpath(@__DIR__, "effective_model", "effective_update.jl"))
include(joinpath(@__DIR__, "effective_model", "effective_evolution.jl"))
include(joinpath(@__DIR__, "effective_model", "effective_horizon_quadratic.jl"))
include(joinpath(@__DIR__, "plots", "boundary_plot.jl"))
include(joinpath(@__DIR__, "plots", "effective_contour.jl"))
include(joinpath(@__DIR__, "plots", "effective_quadratic_horizon_plot.jl"))
include(joinpath(@__DIR__, "plots", "effective_contour_with_quadratic_horizons.jl"))

println("Parameters loaded:")
println("m = ", m)
println("zeta = ", zeta)
println("kappa = ", kappa)
println("alpha = ", alpha)
println("Rc = ", Rc)
println("nCut = ", nCut)
println("z- range = [", zMinusMin, ", ", zMinusMax, "]")
println("z+ range = [", zPlusMin, ", ", zPlusMax, "]")
println("delta z+ = ", deltaZPlus)
println("base z- grid points = ", nMinus)
println("base z- step range = ", extrema(diff(zMinusGrid)))
println("z+ grid points = ", nPlus)
println("outer horizon radius = ", rOuter)
println("outer horizon z- = ", zMinusOuter)
println(
    "automatic z- breakpoints = ",
    (zMinusFarEnd, zMinusMiddleEnd, zMinusFineStart),
)

println()

# Build the zMinus grid used by the effective evolution.
effectiveZMinusGrid = zMinusGrid
effectiveZMinusGrid = denseZMinusGrid(
    effectiveZMinusGrid,
    denseZMinusLower,
    denseZMinusUpper,
    denseDeltaMinus,
)

effectiveZMinusGrid = denseZMinusGrid(
    effectiveZMinusGrid,
    stableZMinusLower,
    stableZMinusUpper,
    stableDeltaMinus,
)
effectiveZMinusGrid = denseZMinusGrid(
    effectiveZMinusGrid,
    ultraStableZMinusLower,
    ultraStableZMinusUpper,
    ultraStableDeltaMinus,
)
effectiveZMinusGrid = denseZMinusGrid(
    effectiveZMinusGrid,
    endpointZMinusLower,
    endpointZMinusUpper,
    endpointDeltaMinus,
)
effectiveZMinusGrid = denseZMinusGrid(
    effectiveZMinusGrid,
    finalZMinusLower,
    finalZMinusUpper,
    finalDeltaMinus,
)
effectiveZMinusGrid = denseZMinusGrid(
    effectiveZMinusGrid,
    lastZMinusLower,
    lastZMinusUpper,
    lastDeltaMinus,
)
effectiveZMinusGrid = denseZMinusGrid(
    effectiveZMinusGrid,
    lastlastZMinusLower,
    lastlastZMinusUpper,
    lastlastDeltaMinus,
)

boundaryData = computeBoundaryData(effectiveZMinusGrid)

using DelimitedFiles
gridOutputFolder = normpath(joinpath(@__DIR__, "output"))

writedlm(
    joinpath(gridOutputFolder, "z_minus_grid.txt"),
    boundaryData.CM.zMinus,
)

writedlm(
    joinpath(gridOutputFolder, "z_plus_grid.txt"),
    boundaryData.CN.zPlus,
)

# cornerData = solveEffectiveCorner(boundaryData = boundaryData, absoluteTolerance = 1.0e-11,
#     relativeTolerance = 1.0e-7)

cornerData =  solveEffectiveCornerAndSaveChunks(boundaryData = boundaryData, saveEveryZPlus = 50.0, absoluteTolerance = 1.0e-11,
    relativeTolerance = 1.0e-7)

cornerDataPath = normpath(joinpath(@__DIR__, "output", "cornerData.jls"))
open(cornerDataPath, "w") do io
    serialize(io, cornerData)
end
println("cornerData saved to: ", cornerDataPath)

horizonData = solveEffectiveQuadraticHorizon(cornerData;zMinusStart = -42.0, halfWindow = 5)
horizonPlot = saveQuadraticHorizonPlot(horizonData; zMinusLimits = (-41.0, -0.1), outputName = "quadratic_horizons.png")

saveEffectiveContourWithQuadraticHorizons(cornerData, horizonData; zMinusLimits = (-41.0, -0.1), zPlusLimits = (0.0, 30.0), levels = 50, outputName = "effective_contours_with_horizons.png");
