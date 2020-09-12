using GeoEstimation
using GeoStatsBase
using Distances
using Plots, VisualRegressionTests
using Test, Pkg, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
isappveyor = "APPVEYOR" ∈ keys(ENV)
isCI = istravis || isappveyor
visualtests = !isCI || (istravis && islinux)
if !isCI
  Pkg.add("Gtk")
  using Gtk
end
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
