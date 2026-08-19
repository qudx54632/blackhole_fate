# Basic functions for the Bardeen vacuum solution.

"""
    RM(zMinus)

Radius on CM: zPlus = 0, normalized by RM(zMinusMin) = Rc.
"""
RM(zMinus) = Rc - (zMinus - zMinusMin) / (2.0 * alpha)

"""
    Phi(r)

Reconstruct the dilaton field Phi = r^2.
"""
Phi(r) = r^2

"""
    Theta(r, J)

Reconstruct Theta = r^2 J from the radial flow and its sensitivity.
"""
Theta(r, J) = r^2 * J

"""
    F(r)

Bardeen radial function.
"""
F(r) = 1.0 - (2.0 * m * r^2) / (r^2 + zeta^2)^(3 / 2)

"""
    Fp(r)

First derivative of F with respect to r.
"""
Fp(r) = (2.0 * m * r * (r^2 - 2.0 * zeta^2)) /
        (r^2 + zeta^2)^(5 / 2)
