# Locate the horizon by fitting a quadratic near the minimum of Phi
# on each fixed-zMinus slice.


# Return the Phi(zPlus) curve on the grid row closest to one chosen zMinus.
function phiCurveAtZMinus(cornerData, zMinusValue)
    rowIndex = argmin(abs.(cornerData.zMinus .- zMinusValue))

    return (
        zMinus = cornerData.zMinus[rowIndex],
        zPlus = cornerData.zPlus,
        phi = cornerData.phiData[rowIndex, :],
    )
end


# Fit a quadratic a*z^2 + b*z + c near the minimum of Phi on one row.
function fitQuadraticPhiMinimum(zPlus, phiValues; halfWindow = 5)
    j0 = argmin(phiValues)
    left = max(1, j0 - halfWindow)
    right = min(length(zPlus), j0 + halfWindow)

    right - left + 1 >= 3 || return nothing

    x = zPlus[left:right]
    y = phiValues[left:right]

    coeffs = hcat(x .^ 2, x, ones(length(x))) \ y
    a, b, c = coeffs

    a > 0.0 || return nothing

    zVertex = -b / (2.0 * a)
    zVertex >= first(x) && zVertex <= last(x) || return nothing

    phiVertex = a * zVertex^2 + b * zVertex + c
    isfinite(zVertex) && isfinite(phiVertex) || return nothing

    return (
        zPlusVertex = zVertex,
        phiVertex = phiVertex,
        a = a,
        b = b,
        c = c,
        minimumIndex = j0,
        window = (left, right),
    )
end


# March from the outer side inward and keep the quadratic minimum on each row.
function solveEffectiveQuadraticHorizon(
    cornerData;
    zMinusStart = -41.0,
    zMinusStop = nothing,
    halfWindow = 5,
)
    startIndex = argmin(abs.(cornerData.zMinus .- zMinusStart))
    stopIndex = zMinusStop === nothing ? length(cornerData.zMinus) :
                                        argmin(abs.(cornerData.zMinus .- zMinusStop))

    if startIndex > stopIndex
        startIndex, stopIndex = stopIndex, startIndex
    end

    points = NamedTuple[]

    for i in startIndex:stopIndex
        fit = fitQuadraticPhiMinimum(
            cornerData.zPlus,
            cornerData.phiData[i, :];
            halfWindow = halfWindow,
        )

        fit === nothing && continue

        push!(points, (; zMinus = cornerData.zMinus[i], fit...))
    end

    return points
end
