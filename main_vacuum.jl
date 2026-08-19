# Main file for the classical vacuum Bardeen solution.

include(joinpath(@__DIR__, "parameters.jl"))
include(joinpath(@__DIR__, "bardeen_vacuum", "functions.jl"))
include(joinpath(@__DIR__, "bardeen_vacuum", "horizons.jl"))
include(joinpath(@__DIR__, "grid.jl"))
include(joinpath(@__DIR__, "bardeen_vacuum", "radial_flow.jl"))
include(joinpath(@__DIR__, "plots", "vacuum_contour.jl"))

println("Calculating the classical vacuum Bardeen solution...")
vacuumData = solveVacuumRegion()
classicalHorizons = calculateHorizons()

if classicalHorizons === nothing
    println("No classical Bardeen horizons exist for these parameters.")
else
    println("Inner horizon:")
    println("    r = ", classicalHorizons.inner.radius)
    println("    zMinus = ", classicalHorizons.inner.zMinus)

    println("Outer horizon:")
    println("    r = ", classicalHorizons.outer.radius)
    println("    zMinus = ", classicalHorizons.outer.zMinus)
end

vacuumPlotPath = saveVacuumContour(vacuumData, classicalHorizons)
println("Vacuum contour plot saved to: ", vacuumPlotPath)
