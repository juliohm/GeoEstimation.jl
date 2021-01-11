# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    IDW(var₁=>param₁, var₂=>param₂, ...)

Inverse distance weighting estimation solver.

## Parameters

* `neighbors` - Number of neighbors (default to all the data)
* `distance`  - A distance defined in Distances.jl (default to `Euclidean()`)
* `power`     - Power of the distances (default to `1`)

### References

Shepard 1968. *A two-dimensional interpolation function for irregularly-spaced data.*
"""
@estimsolver IDW begin
  @param neighbors = nothing
  @param distance = Euclidean()
  @param power = 1
end

function solve(problem::EstimationProblem, solver::IDW)
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
      n = length(z)

      # determine number of nearest neighbors to use
      k = isnothing(varparams.neighbors) ? n : varparams.neighbors

      # determine distance type
      D = varparams.distance

      # determine power of distances
      p = varparams.power

      @assert n > 0 "estimation requires data"

      @assert k ≤ n "invalid number of neighbors"

      @assert p > 0 "power must be positive"

      # fit search tree
      if D isa NearestNeighbors.MinkowskiMetric
        tree = KDTree(X, D)
      else
        tree = BallTree(X, D)
      end

      # pre-allocate memory for results
      varμ = Vector{V}(undef, nelms(pdomain))
      varσ = Vector{V}(undef, nelms(pdomain))

      # pre-allocate memory for coordinates
      x = MVector{N,T}(undef)

      # estimation loop
      for loc in traverse(pdomain, LinearPath())
        coordinates!(x, pdomain, loc)

        is, ds = knn(tree, x, k)
        ws = one(V) ./ ds.^p
        Σw = sum(ws)

        if isinf(Σw) # some distance is zero?
          j = findfirst(iszero, ds)
          μ = z[is[j]]
          σ = zero(V)
        else
          ws /= Σw
          vs  = view(z, is)
          μ = sum(ws[i]*vs[i] for i in eachindex(vs))
          σ = minimum(ds)
        end

        varμ[loc] = μ
        varσ[loc] = σ
      end

      push!(μs, var => varμ)
      push!(σs, Symbol(var,:Var) => varσ)
    end
  end

  georef((; μs..., σs...), pdomain)
end
