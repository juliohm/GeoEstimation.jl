using GeoEstimation
using Meshes
using GeoStatsBase
using Variography
using Distances
using Plots; gr(size=(600,400))
using ReferenceTests, ImageIO
using Test, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__,"data")

# list of tests
testfiles = [
  "idw.jl",
  "lwr.jl",
  "kriging.jl"
]

@testset "GeoEstimation.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
