# Bardeen horizon radii and their double-null coordinate locations.

using Roots

"""
    horizonZMinus(radius)

Return the value of zMinus where RM(zMinus) = radius.
"""
horizonZMinus(radius) =
    zMinusMin + 2.0 * alpha * (Rc - radius)

"""
    calculateHorizons()

Find the inner and outer positive roots of F(r) = 0 and their
constant-zMinus locations. Return nothing when no horizon exists.
"""
function calculateHorizons()
    rCritical = sqrt(2.0) * zeta
    minimumF = F(rCritical)

    if minimumF > 0.0
        return nothing
    end

    if isapprox(minimumF, 0.0; atol = 1.0e-12, rtol = 1.0e-12)
        horizon = (radius = rCritical, zMinus = horizonZMinus(rCritical))
        return (inner = horizon, outer = horizon)
    end

    rInner = find_zero(F, (0.0, rCritical), Bisection())

    rRight = max(2.0 * m, 2.0 * rCritical)
    while F(rRight) <= 0.0
        rRight *= 2.0
    end
    rOuter = find_zero(F, (rCritical, rRight), Bisection())

    return (
        inner = (radius = rInner, zMinus = horizonZMinus(rInner)),
        outer = (radius = rOuter, zMinus = horizonZMinus(rOuter)),
    )
end
