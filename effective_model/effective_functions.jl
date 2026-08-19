# Local effective-model formulas.
#
# Stored variables:
#
#     Z = log(Phi),   X = log(Theta / Phi).


# ---------------------------------------------------------------------------
# Variable conversions
# ---------------------------------------------------------------------------

Z(phi) = log(phi)
X(theta, phi) = log(theta / phi)

phiFromZ(z) = exp(z)
thetaFromZX(z, x) = exp(z + x)


# ---------------------------------------------------------------------------
# Basic positivity check
# ---------------------------------------------------------------------------

function requirePositivePhi(phi)
    phi > 0.0 ||
        throw(DomainError(phi, "Phi must be positive in the effective equations."))

    return nothing
end


# ---------------------------------------------------------------------------
# Effective-equation coefficients
# ---------------------------------------------------------------------------

# Coefficient multiplying P = ∂+∂-Phi.
function A(phi)
    requirePositivePhi(phi)
    return (phi + zeta^2)^(3 / 2) / (4.0 * phi^(3 / 2))
end

# Constant coefficient multiplying Y = ∂+∂-X.
b() = kappa * zeta^2 / 24.0

# Phi-dependent coefficient multiplying Y.
function c(phi)
    requirePositivePhi(phi)
    return (
        kappa * zeta^4 * sqrt(phi) /
        (2.0 * (phi + zeta^2)^(5 / 2))
    )
end


# ---------------------------------------------------------------------------
# Source terms and determinant
# ---------------------------------------------------------------------------

# Source term in the first effective equation.
function B(phi, x, pMinusValue, pPlusValue)
    requirePositivePhi(phi)

    zetaSquared = zeta^2
    derivativeProduct = pMinusValue * pPlusValue

    return (
        phi * (phi - 2.0 * zetaSquared) * exp(x) -
        3.0 * zetaSquared * derivativeProduct
    ) / (2.0 * phi * (phi + zetaSquared))
end

# Source term in the second effective equation.
function STheta(phi, x, pMinusValue, pPlusValue)
    requirePositivePhi(phi)

    zetaSquared = zeta^2
    derivativeProduct = pMinusValue * pPlusValue

    return (
        (
            13.0 * zetaSquared * phi +
            2.0 * phi^2 -
            4.0 * zetaSquared^2
        ) * derivativeProduct +
        3.0 * zetaSquared * phi * exp(x) *
        (3.0 * phi - 2.0 * zetaSquared)
    ) / (4.0 * phi^2 * (phi + zetaSquared)^2)
end

# Reduced determinant used in the stable formula for Y.
function DeltaTilde(phi)
    requirePositivePhi(phi)
    return phi * (1.0 + c(phi)) - b() / A(phi)
end

# Geometry factor used in the stable formula for Y.
function G(phi)
    requirePositivePhi(phi)

    zetaSquared = zeta^2
    return 2.0 * phi^2 - 11.0 * zetaSquared * phi + 2.0 * zetaSquared^2
end


# ---------------------------------------------------------------------------
# Solve the local algebraic system for P and Y
# ---------------------------------------------------------------------------

"""
    solvePY(phi, x, pMinus, pPlus)

Solve the two local effective equations at one diamond center.
"""
function solvePY(phi, x, pMinusValue, pPlusValue)
    requirePositivePhi(phi)

    aCoefficient = A(phi)
    bCoefficient = b()
    deltaTilde = DeltaTilde(phi)

    regularCenterCombination = pMinusValue * pPlusValue + phi * exp(x)

    Y =
        G(phi) * regularCenterCombination /
        (
            4.0 * phi * (phi + zeta^2)^2 *
            deltaTilde
        )

    P = -B(phi, x, pMinusValue, pPlusValue) -
        (bCoefficient / aCoefficient) * Y

    isfinite(P) && isfinite(Y) ||
        error("P or Y is not finite.")

    return (P = P, Y = Y, DeltaTilde = deltaTilde)
end
