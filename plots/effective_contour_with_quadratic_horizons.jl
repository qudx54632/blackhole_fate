# Effective contour plots with quadratic-fit horizons overlaid.

using Plots
using LaTeXStrings


function addQuadraticHorizonCurves!(fieldPlot, branches)
    plot!(
        fieldPlot,
        branches.outerHorizon.zPlus,
        branches.outerHorizon.zMinus;
        color = :red,
        linewidth = 2.5,
        label = L"\mathrm{outer\ horizon}",
    )

    plot!(
        fieldPlot,
        branches.innerHorizon.zPlus,
        branches.innerHorizon.zMinus;
        color = :green,
        linewidth = 2.5,
        label = L"\mathrm{inner\ horizon}",
    )

    return fieldPlot
end


function saveEffectiveContourWithQuadraticHorizons(
    cornerData,
    horizonData;
    zMinusLimits = extrema(cornerData.zMinus),
    zPlusLimits = extrema(cornerData.zPlus),
    levels = 50,
    outputName = "effective_contours_with_quadratic_horizons.png",
    splitIndex = nothing,
)
    minusIndices = findall(
        zMinus -> zMinusLimits[1] <= zMinus <= zMinusLimits[2],
        cornerData.zMinus,
    )
    plusIndices = findall(
        zPlus -> zPlusLimits[1] <= zPlus <= zPlusLimits[2],
        cornerData.zPlus,
    )

    zMinus = cornerData.zMinus[minusIndices]
    zPlus = cornerData.zPlus[plusIndices]
    phiData = cornerData.phiData[minusIndices, plusIndices]
    thetaData = cornerData.thetaData[minusIndices, plusIndices]

    branches = splitQuadraticHorizonBranches(
        horizonData;
        splitIndex = splitIndex,
    )

    phiPlot = contour(
        zPlus,
        zMinus,
        phiData;
        fill = true,
        levels = levels,
        color = :viridis,
        xlabel = L"z^+",
        ylabel = L"z^-",
        title = L"\Phi_{\mathrm{effective}}",
        xlims = zPlusLimits,
        ylims = zMinusLimits,
        colorbar = true,
        legend = :best,
    )
    addQuadraticHorizonCurves!(phiPlot, branches)

    thetaPlot = contour(
        zPlus,
        zMinus,
        thetaData;
        fill = true,
        levels = levels,
        color = :viridis,
        xlabel = L"z^+",
        ylabel = L"z^-",
        title = L"\Theta_{\mathrm{effective}}",
        xlims = zPlusLimits,
        ylims = zMinusLimits,
        colorbar = true,
        legend = :best,
    )
    addQuadraticHorizonCurves!(thetaPlot, branches)

    figure = plot(
        phiPlot,
        thetaPlot;
        layout = (1, 2),
        size = (1400, 550),
        left_margin = 5Plots.mm,
        bottom_margin = 6Plots.mm,
    )

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
