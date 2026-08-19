# Optional full-vacuum calculation used to inspect the boundary construction.

using Plots
using LaTeXStrings

"""
    solveVacuumRegion()

Solve the radial-flow ODEs for every value of zMinus and reconstruct
Phi and Theta on the full vacuum double-null grid.
"""
function solveVacuumRegion()
    zMinus = copy(zMinusGrid)
    zPlus = collect(range(zPlusMin, zPlusMax, length = nPlus))

    rData = Matrix{Float64}(undef, nPlus, nMinus)
    JData = Matrix{Float64}(undef, nPlus, nMinus)

    for (i, zm) in enumerate(zMinus)
        solution = solveRadialFlow(RM(zm))
        rData[:, i] .= [solution(zp; idxs = 1) for zp in zPlus]
        JData[:, i] .= [solution(zp; idxs = 2) for zp in zPlus]
    end

    phiData = Phi.(rData)
    thetaData = Theta.(rData, JData)

    return (; zMinus, zPlus, rData, JData, phiData, thetaData)
end

"""
    saveVacuumContour(data, horizons)

Save side-by-side vacuum contour plots of Phi and Theta, with the inner
and outer horizons overlaid when they exist.
"""
function saveVacuumContour(data, horizons)
    # plotIndices = findall(
    #     zMinus -> zMinusFineStart <= zMinus <= zMinusMax,
    #     data.zMinus,
    # )

    plotIndices = findall(
        zMinus -> zMinusFineStart <= zMinus <= -35.0,
        data.zMinus,
    )

    plotZMinus = data.zMinus[plotIndices]
    plotPhi = data.phiData[:, plotIndices]
    plotTheta = data.thetaData[:, plotIndices]

    phiPlot = contour(
        data.zPlus,
        plotZMinus,
        permutedims(plotPhi);
        fill = true,
        levels = 100,
        color = :viridis,
        xlabel = L"z^+",
        ylabel = L"z^-",
        title = L"\Phi",
        colorbar = true,
        legend = :topright,
    )

    thetaPlot = contour(
        data.zPlus,
        plotZMinus,
        permutedims(plotTheta);
        fill = true,
        levels = 100,
        color = :viridis,
        xlabel = L"z^+",
        ylabel = L"z^-",
        title = L"\Theta",
        colorbar = true,
        legend = :topright,
    )

    if horizons !== nothing
        for fieldPlot in (phiPlot, thetaPlot)
            # hline!(
            #     fieldPlot,
            #     [horizons.inner.zMinus];
            #     color = :orange,
            #     linestyle = :dashdot,
            #     linewidth = 2,
            #     label = L"\mathrm{Inner\ horizon}",
            # )
            hline!(
                fieldPlot,
                [horizons.outer.zMinus];
                color = :red,
                linestyle = :dash,
                linewidth = 2,
                label = L"\mathrm{Outer\ horizon}",
            )
        end
    end

    figure = plot(phiPlot, thetaPlot; layout = (1, 2), size = (1400, 550))
    outputPath = normpath(joinpath(@__DIR__, "..", "output", "vacuum_contours.png"))
    mkpath(dirname(outputPath))
    savefig(figure, outputPath)

    return outputPath
end
