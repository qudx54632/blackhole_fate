# Boundary values on the two initial null lines.

function computeBoundaryData(zMinusValues = zMinusGrid)
    zMinus = Float64.(zMinusValues)
    zPlus = collect(range(zPlusMin, zPlusMax, length = nPlus))

    # CM: zPlus = 0.
    rM = RM.(zMinus)
    phiM = rM .^ 2

    CM = (
        zMinus = zMinus,
        R = rM,
        Phi = phiM,
        Theta = phiM,
        X = zeros(length(zMinus)),
    )

    # CN: zMinus = zMinusMin.
    soln = solveRadialFlow(Rc)
    rN = [soln(zp; idxs = 1) for zp in zPlus]
    JN = [soln(zp; idxs = 2) for zp in zPlus]
    phiN = rN .^ 2

    @assert all(JN .> 0.0) "J must be positive because X = log(J)."

    CN = (
        zPlus = zPlus,
        R = rN,
        J = JN,
        Phi = phiN,
        Theta = phiN .* JN,
        X = log.(JN),
    )

    return (CM = CM, CN = CN)
end
