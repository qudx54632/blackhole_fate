# Plot the quadratic-fit horizons.

using Plots
using LaTeXStrings


function quadraticHorizonCoordinates(horizonData)
    zMinus = [horizonData[i].zMinus for i in 1:length(horizonData)]
    zPlus = [horizonData[i].zPlusVertex for i in 1:length(horizonData)]
    phi = [horizonData[i].phiVertex for i in 1:length(horizonData)]

    return zMinus, zPlus, phi
end


function orderHorizonBranch(zMinus, zPlus)
    order = sortperm(zPlus)
    return (
        zMinus = zMinus[order],
        zPlus = zPlus[order],
    )
end


function splitQuadraticHorizonBranches(horizonData; splitIndex = nothing)
    zMinus, zPlus, phi = quadraticHorizonCoordinates(horizonData)

    order = sortperm(zMinus)
    zMinus = zMinus[order]
    zPlus = zPlus[order]
    phi = phi[order]

    if splitIndex === nothing
        gaps = abs.(diff(phi))
        splitIndex = argmax(gaps)
    end
    splitPhiJump = abs(phi[splitIndex + 1] - phi[splitIndex])

    zMinusA = zMinus[1:splitIndex]
    zPlusA = zPlus[1:splitIndex]
    zMinusB = zMinus[(splitIndex + 1):end]
    zPlusB = zPlus[(splitIndex + 1):end]

    meanMinusA = sum(zMinusA) / length(zMinusA)
    meanMinusB = sum(zMinusB) / length(zMinusB)

    if meanMinusA > meanMinusB
        innerHorizon = orderHorizonBranch(zMinusA, zPlusA)
        outerHorizon = orderHorizonBranch(zMinusB, zPlusB)
    else
        innerHorizon = orderHorizonBranch(zMinusB, zPlusB)
        outerHorizon = orderHorizonBranch(zMinusA, zPlusA)
    end

    return (
        innerHorizon = innerHorizon,
        outerHorizon = outerHorizon,
        splitIndex = splitIndex,
        splitPhiJump = splitPhiJump,
    )
end


function saveQuadraticHorizonPlot(
    horizonData;
    outputName = "effective_quadratic_horizons.png",
    zMinusLimits = nothing,
    zPlusLimits = nothing,
    splitIndex = nothing,
)
    branches = splitQuadraticHorizonBranches(
        horizonData;
        splitIndex = splitIndex,
    )

    figure = plot(
        branches.outerHorizon.zPlus,
        branches.outerHorizon.zMinus;
        color = :red,
        linewidth = 2.5,
        label = L"\mathrm{outer\ horizon}",
        xlabel = L"z^+",
        ylabel = L"z^-",
        title = L"\mathrm{Quadratic\ horizon\ finder}",
        legend = :best,
        grid = true,
    )

    plot!(
        figure,
        branches.innerHorizon.zPlus,
        branches.innerHorizon.zMinus;
        color = :green,
        linewidth = 2.5,
        label = L"\mathrm{inner\ horizon}",
    )

    zMinusLimits === nothing || plot!(figure; ylims = zMinusLimits)
    zPlusLimits === nothing || plot!(figure; xlims = zPlusLimits)

    outputPath = normpath(joinpath(@__DIR__, "..", "output", outputName))
    mkpath(dirname(outputPath))
    savefig(figure, outputPath)

    return (
        innerHorizon = branches.innerHorizon,
        outerHorizon = branches.outerHorizon,
        splitIndex = branches.splitIndex,
        splitPhiJump = branches.splitPhiJump,
        plotPath = outputPath,
    )
end
