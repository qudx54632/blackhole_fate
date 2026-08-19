# Evolve the effective equations with a fixed zPlus step.
# The zMinus grid is fixed during the whole evolution.

using Serialization


# Calculate Z = log(Phi) and X on one new fixed-zPlus slice.
function evolveOneZPlusColumn(
    ZOld,
    xOld,
    zMinus,
    zPlusOld,
    zPlusNew,
    ZBoundary,
    xBoundary;
    absoluteTolerance = 1.0e-11,
    relativeTolerance = 1.0e-7,
)
    nMinus = length(zMinus)
    ZNew = fill(NaN, nMinus)
    xNew = fill(NaN, nMinus)
    deltaPlus = zPlusNew - zPlusOld

    ZNew[1] = ZBoundary
    xNew[1] = xBoundary

    for i in 1:(nMinus - 1)
        deltaMinus = zMinus[i + 1] - zMinus[i]

        south = [ZOld[i], xOld[i]]
        west = [ZOld[i + 1], xOld[i + 1]]
        east = [ZNew[i], xNew[i]]

        zMinusCenter = (zMinus[i] + zMinus[i + 1]) / 2
        zPlusCenter = zPlusOld + deltaPlus / 2

        step = effectiveDiamondStep(
            south,
            east,
            west,
            zMinusCenter,
            zPlusCenter,
            deltaMinus,
            deltaPlus;
            absoluteTolerance = absoluteTolerance,
            relativeTolerance = relativeTolerance,
            maxIterations = 100,
        )

        ZNew[i + 1] = step.Z
        xNew[i + 1] = step.x
    end

    return ZNew, xNew
end


function solveEffectiveCorner(
    ;
    boundaryData = computeBoundaryData(),
    absoluteTolerance = 1.0e-11,
    relativeTolerance = 1.0e-7,
)
    zMinus = boundaryData.CM.zMinus

    # Initial slice on CM: zPlus = 0.
    phiInitial = boundaryData.CM.Phi
    thetaInitial = boundaryData.CM.Theta
    ZColumns = [Z.(phiInitial)]
    xColumns = [X.(thetaInitial, phiInitial)]

    # Boundary values on CN: zMinus = zMinusMin.
    zPlusAll = boundaryData.CN.zPlus
    phiBoundary = boundaryData.CN.Phi
    thetaBoundary = boundaryData.CN.Theta
    ZBoundary = Z.(phiBoundary)
    xBoundary = X.(thetaBoundary, phiBoundary)

    zPlus = Float64[zPlusAll[1]]

    for j in 1:(length(zPlusAll) - 1)
        zPlusOld = zPlusAll[j]
        zPlusNew = zPlusAll[j + 1]

        ZNew, xNew = evolveOneZPlusColumn(
            ZColumns[end],
            xColumns[end],
            zMinus,
            zPlusOld,
            zPlusNew,
            ZBoundary[j + 1],
            xBoundary[j + 1];
            absoluteTolerance = absoluteTolerance,
            relativeTolerance = relativeTolerance,
        )

        push!(ZColumns, ZNew)
        push!(xColumns, xNew)
        push!(zPlus, zPlusNew)
    end

    ZData = hcat(ZColumns...)
    xData = hcat(xColumns...)
    phiData = phiFromZ.(ZData)
    thetaData = thetaFromZX.(ZData, xData)

    return (; zMinus, zPlus, ZData, phiData, xData, thetaData)
end


function saveEffectiveCornerChunk(
    chunkFolder,
    chunkNumber,
    zMinus,
    zPlusSegment,
    ZColumns,
    xColumns,
)
    ZData = hcat(ZColumns...)
    xData = hcat(xColumns...)
    phiData = phiFromZ.(ZData)
    thetaData = thetaFromZX.(ZData, xData)

    chunkName = "corner_chunk_" * lpad(string(chunkNumber), 4, '0') * ".jls"
    chunkPath = joinpath(chunkFolder, chunkName)

    open(chunkPath, "w") do io
        serialize(
            io,
            (
                zMinus = collect(zMinus),
                zPlus = collect(zPlusSegment),
                ZData = ZData,
                phiData = phiData,
                xData = xData,
                thetaData = thetaData,
            ),
        )
    end

    return chunkPath
end


function solveEffectiveCornerAndSaveChunks(
    ;
    boundaryData = computeBoundaryData(),
    saveEveryZPlus = 50.0,
    chunkFolder = normpath(
        joinpath(@__DIR__, "..", "output", "effective_checkpoints"),
    ),
    absoluteTolerance = 1.0e-11,
    relativeTolerance = 1.0e-7,
)
    zMinus = boundaryData.CM.zMinus

    phiInitial = boundaryData.CM.Phi
    thetaInitial = boundaryData.CM.Theta
    ZColumns = [Z.(phiInitial)]
    xColumns = [X.(thetaInitial, phiInitial)]

    zPlusAll = boundaryData.CN.zPlus
    phiBoundary = boundaryData.CN.Phi
    thetaBoundary = boundaryData.CN.Theta
    ZBoundary = Z.(phiBoundary)
    xBoundary = X.(thetaBoundary, phiBoundary)

    mkpath(chunkFolder)
    for oldFile in filter(
        file -> startswith(basename(file), "corner_chunk_") &&
                endswith(file, ".jls"),
        readdir(chunkFolder; join = true),
    )
        rm(oldFile)
    end

    chunkPaths = String[]
    chunkNumber = 0
    segmentZPlus = Float64[zPlusAll[1]]
    nextSaveZPlus = zPlusAll[1] + saveEveryZPlus

    for j in 1:(length(zPlusAll) - 1)
        zPlusOld = zPlusAll[j]
        zPlusNew = zPlusAll[j + 1]

        ZNew, xNew = evolveOneZPlusColumn(
            ZColumns[end],
            xColumns[end],
            zMinus,
            zPlusOld,
            zPlusNew,
            ZBoundary[j + 1],
            xBoundary[j + 1];
            absoluteTolerance = absoluteTolerance,
            relativeTolerance = relativeTolerance,
        )

        push!(ZColumns, ZNew)
        push!(xColumns, xNew)
        push!(segmentZPlus, zPlusNew)

        if zPlusNew >= nextSaveZPlus
            chunkNumber += 1
            push!(
                chunkPaths,
                saveEffectiveCornerChunk(
                    chunkFolder,
                    chunkNumber,
                    zMinus,
                    segmentZPlus,
                    ZColumns,
                    xColumns,
                ),
            )

            ZColumns = [ZColumns[end]]
            xColumns = [xColumns[end]]
            segmentZPlus = Float64[segmentZPlus[end]]

            while zPlusNew >= nextSaveZPlus
                nextSaveZPlus += saveEveryZPlus
            end
        end
    end

    if length(segmentZPlus) > 1 || isempty(chunkPaths)
        chunkNumber += 1
        push!(
            chunkPaths,
            saveEffectiveCornerChunk(
                chunkFolder,
                chunkNumber,
                zMinus,
                segmentZPlus,
                ZColumns,
                xColumns,
            ),
        )
    end

    return (; zMinus, zPlus = zPlusAll, chunkPaths)
end
