# Physical parameters
m = 10.0
zeta = 1.0

# Effective-model quantum coupling. Adjust this value as needed.
kappa = 50.0

# Double-null coordinate normalization
alpha = 1.0

# Corner radius and coordinate cutoff
Rc = 1000.0
nCut = 2.0 * alpha * Rc

# Computational domain
zMinusMin = -nCut
zMinusMax = -0.1
zPlusMin = 0.0
zPlusMax = 200.0
deltaZPlus = 0.01

# Number of points on the boundary zMinus = zMinusMin.
nPlus = round(Int, (zPlusMax - zPlusMin) / deltaZPlus) + 1

# Each zMinus breakpoint is placed where RM is the specified multiple of
# the automatically calculated outer-horizon radius.
gridFarRadiusFactor = 5.0
gridMiddleRadiusFactor = 1.5
gridFineRadiusFactor = 1.125

# Resolution in each part of the piecewise zMinus grid.
deltaMinusFar = 10.0
deltaMinusMiddle = 1.0
deltaMinusTransition = 0.01
deltaMinusFine = 0.01

# # For the effective evolution, use a moderately fine zMinus corridor around
# # the horizon region, and a very fine core only near the late-time stable part.
denseZMinusLower = -39.20
denseZMinusUpper = -36.80
denseDeltaMinus = 0.005

stableZMinusLower = -39.05
stableZMinusUpper = -38.00
stableDeltaMinus = 0.0005

# If the late-time horizon becomes very stiff, make one much narrower strip
# even finer.  The present failures are around zMinus = -38.89.
ultraStableZMinusLower = -38.98
ultraStableZMinusUpper = -38.86
ultraStableDeltaMinus = 0.0001

# Very narrow endpoint core near the present late-time failure point.
endpointZMinusLower = -38.91
endpointZMinusUpper = -38.90
endpointDeltaMinus = 1.0e-5

# Extremely narrow core for the later stiff points near zMinus = -38.909.
finalZMinusLower = -38.9092
finalZMinusUpper = -38.9080
finalDeltaMinus = 1.0e-6

# Last very narrow core near the late-time endpoint, currently around
# zMinus = -38.9088945.
lastZMinusLower = -38.90898
lastZMinusUpper = -38.90884
lastDeltaMinus = 1.0e-8

lastlastZMinusLower = -38.90898
lastlastZMinusUpper = -38.90893
lastlastDeltaMinus = 1.0e-8


