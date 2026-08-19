# Plot Phi and Theta on the two characteristic boundaries.

using Plots
using LaTeXStrings

"""
    saveBoundaryDataPlot(boundaryData)

Save one overlaid Phi/Theta panel for each of CM and CN.
"""
function saveBoundaryDataPlot(boundaryData)
    boundaryMPlot = plot(
        boundaryData.CM.zMinus,
        boundaryData.CM.Phi;
        color = :navy,
        linewidth = 3,
        label = L"\Phi",
        xlabel = L"z^-",
        ylabel = L"\mathrm{Field\ value}",
        title = L"\mathcal{C}_M:\ z^+=0",
        xflip = true,
        legend = :topleft,
        grid = true,
    )
    plot!(
        boundaryMPlot,
        boundaryData.CM.zMinus,
        boundaryData.CM.Theta;
        color = :darkorange,
        linewidth = 2,
        linestyle = :dash,
        label = L"\Theta",
    )

    boundaryNPlot = plot(
        boundaryData.CN.zPlus,
        boundaryData.CN.Phi;
        color = :navy,
        linewidth = 3,
        label = L"\Phi",
        xlabel = L"z^+",
        ylabel = L"\mathrm{Field\ value}",
        title = L"\mathcal{C}_N:\ z^-=-N",
        legend = :topleft,
        grid = true,
    )
    plot!(
        boundaryNPlot,
        boundaryData.CN.zPlus,
        boundaryData.CN.Theta;
        color = :darkorange,
        linewidth = 2,
        linestyle = :dash,
        label = L"\Theta",
    )

    figure = plot(
        boundaryMPlot,
        boundaryNPlot;
        layout = (1, 2),
        size = (1200, 550),
        bottom_margin = 8Plots.mm,
        left_margin = 6Plots.mm,
        plot_title = L"\mathrm{Characteristic\ boundary\ data}",
    )

    outputPath = normpath(
        joinpath(@__DIR__, "..", "output", "boundary_data.png"),
    )
    mkpath(dirname(outputPath))
    savefig(figure, outputPath)

    return outputPath
end
