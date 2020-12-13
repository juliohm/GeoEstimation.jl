# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    LWR(var₁=>param₁, var₂=>param₂, ...)

Locally weighted regression estimation solver.

## Parameters

* `weightfun` - Weighting function (default to `exp(-h^2/2)`)
* `distance`  - A distance from Distances.jl (default to `Euclidean()`)
* `neighbors` - Number of neighbors (default to 20% of the data)

### References

* Stone 1977. *Consistent non-parametric regression.*
* Cleveland 1979. *Robust locally weighted regression and smoothing scatterplots.*
* Cleveland & Grosse 1991. *Computational methods for local regression.*
"""
@estimsolver LWR begin
  @param weightfun = h -> exp(-3*h^2)
  @param distance = Euclidean()
  @param neighbors = nothing
end

function solve(problem::EstimationProblem, solver::LWR)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)
  N = ncoords(pdomain)
  T = coordtype(pdomain)

  mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

  # result for each variable
  μs = []; σs = []

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine value type
      V = mactypeof[var]

      # retrieve non-missing data
      locs = findall(!ismissing, pdata[var])
      𝒟 = view(pdata, locs)
      X = coordinates(𝒟)
      z = 𝒟[var]

      # number of data points for variable
      ndata = length(z)

      # weight function
      w = varparams.weightfun

      # number of nearest neighbors
      k = isnothing(varparams.neighbors) ? ceil(Int, 0.2ndata) : varparams.neighbors

      @assert 0 < k ≤ ndata "invalid number of neighbors"

      # fit search tree
      M = varparams.distance
      if M isa NearestNeighbors.MinkowskiMetric
        tree = KDTree(X, M)
      else
        tree = BallTree(X, M)
      end

      # pre-allocate memory for results
      varμ = Vector{V}(undef, nelms(pdomain))
      varσ = Vector{V}(undef, nelms(pdomain))

      # pre-allocate memory for coordinates
      x = MVector{N,T}(undef)

      # estimation loop
      for loc in traverse(pdomain, LinearPath())
        coordinates!(x, pdomain, loc)

        # find neighbors
        is, ds = knn(tree, x, k)
        δs = ds ./ maximum(ds)

        # weighted least-squares
        Wₗ = Diagonal(w.(δs))
        Xₗ = [ones(eltype(X), k) X[:,is]']
        zₗ = view(z, is)
        θₗ = Xₗ'*Wₗ*Xₗ \ Xₗ'*Wₗ*zₗ

        # linear combination of response values
        xₒ = [one(eltype(x)); x]
        ẑₒ = θₗ ⋅ xₒ
        rₗ = Wₗ*Xₗ*(Xₗ'*Wₗ*Xₗ\xₒ)
        r̂ₒ = norm(rₗ)

        varμ[loc] = ẑₒ
        varσ[loc] = r̂ₒ
      end

      push!(μs, var => varμ)
      push!(σs, var => varσ)
    end
  end

  EstimationSolution(pdomain, Dict(μs), Dict(σs))
end
