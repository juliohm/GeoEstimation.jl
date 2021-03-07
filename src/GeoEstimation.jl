# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

module GeoEstimation

using Meshes
using GeoStatsBase
using Variography
using KrigingEstimators

using NearestNeighbors
using Distances
using LinearAlgebra

import GeoStatsBase: solve

include("idw.jl")
include("lwr.jl")
include("kriging.jl")

export IDW, LWR, Kriging

end
