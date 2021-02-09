using GeoEstimation
using GeoStatsBase
using Distances
using Plots, VisualRegressionTests
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
  "lwr.jl"
]

@testset "GeoEstimation.jl" begin
  for testfile in testfiles
    include(testfile)
  end
end
