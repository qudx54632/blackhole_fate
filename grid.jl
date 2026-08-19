# Construct the fixed, nonuniform zMinus grid.

# Classical outer horizon, used only to place the initial grid.
classicalHorizonsForGrid = calculateHorizons()


rOuter = classicalHorizonsForGrid.outer.radius
zMinusOuter = classicalHorizonsForGrid.outer.zMinus

# Radii at which the grid resolution changes.
rFar = gridFarRadiusFactor * rOuter
rMiddle = gridMiddleRadiusFactor * rOuter
rFine = gridFineRadiusFactor * rOuter

# The same three points in the zMinus coordinate.
zMinusFarEnd = horizonZMinus(rFar)
zMinusMiddleEnd = horizonZMinus(rMiddle)
zMinusFineStart = horizonZMinus(rFine)

# Number of cells in each region. Ceil ensures that the actual step is no
# larger than the requested deltaMinus value.
nFarCells = ceil(Int, (zMinusFarEnd - zMinusMin) / deltaMinusFar)
nMiddleCells = ceil(
    Int,
    (zMinusMiddleEnd - zMinusFarEnd) / deltaMinusMiddle,
)
nTransitionCells = ceil(
    Int,
    (zMinusFineStart - zMinusMiddleEnd) / deltaMinusTransition,
)
nFineCells = ceil(Int, (zMinusMax - zMinusFineStart) / deltaMinusFine)

# Uniform points inside each of the four regions.
zMinusFarGrid = collect(
    range(zMinusMin, zMinusFarEnd, length = nFarCells + 1),
)
zMinusMiddleGrid = collect(
    range(zMinusFarEnd, zMinusMiddleEnd, length = nMiddleCells + 1),
)
zMinusTransitionGrid = collect(
    range(
        zMinusMiddleEnd,
        zMinusFineStart,
        length = nTransitionCells + 1,
    ),
)
zMinusFineGrid = collect(
    range(zMinusFineStart, zMinusMax, length = nFineCells + 1),
)

# Join the four regions. The last point of one region is the first point of
# the next, so [1:end-1] avoids duplicate points at the interfaces.
zMinusGrid = vcat(
    zMinusFarGrid[1:end-1],
    zMinusMiddleGrid[1:end-1],
    zMinusTransitionGrid[1:end-1],
    zMinusFineGrid,
)

nMinus = length(zMinusGrid)


# Divide each zMinus interval inside one chosen band into smaller intervals.
function refineZMinusGrid(zMinus, lower, upper, factor)
    refinedZMinus = Float64[first(zMinus)]

    for i in 1:(length(zMinus) - 1)
        left = zMinus[i]
        right = zMinus[i + 1]

        if right >= lower && left <= upper
            newPoints = range(left, right, length = factor + 1)
            append!(refinedZMinus, newPoints[2:end])
        else
            push!(refinedZMinus, right)
        end
    end

    return refinedZMinus
end


# Replace one zMinus band by a uniformly dense grid.
function denseZMinusGrid(zMinus, lower, upper, delta)
    lower < upper || error("The dense zMinus interval must have lower < upper.")
    delta > 0.0 || error("denseDeltaMinus must be positive.")

    leftGrid = zMinus[zMinus .< lower]
    rightGrid = zMinus[zMinus .> upper]

    numberOfDenseCells = ceil(Int, (upper - lower) / delta)
    denseGrid = collect(
        range(lower, upper, length = numberOfDenseCells + 1),
    )

    return sort(unique(vcat(leftGrid, denseGrid, rightGrid)))
end
