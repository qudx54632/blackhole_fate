# Contour plots for the effective-model fields.

using Plots
using LaTeXStrings


function saveEffectiveContourPlots(
    cornerData;
    zMinusLimits = extrema(cornerData.zMinus),
    zPlusLimits = extrema(cornerData.zPlus),
    levels = 50,
    outputName = "effective_contours.png",
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
        colorbar = true,
    )

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
        colorbar = true,
    )

    figure = plot(
        phiPlot,
        thetaPlot;
        layout = (1, 2),
        size = (1300, 550),
        left_margin = 5Plots.mm,
        bottom_margin = 6Plots.mm,
    )

    outputPath = normpath(joinpath(@__DIR__, "..", "output", outputName))
    mkpath(dirname(outputPath))
    savefig(figure, outputPath)

    return outputPath
end
