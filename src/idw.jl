# ------------------------------------------------------------------
# Licensed under the ISC License. See LICENCE in the project root.
# ------------------------------------------------------------------

"""
    InvDistWeight(var₁=>param₁, var₂=>param₂, ...)

Inverse distance weighting estimation solver.

## Parameters

* `neighbors` - Number of neighbors (default to all data locations)
* `distance`  - A distance defined in Distances.jl (default to `Euclidean()`)

### References

Shepard 1968. *A two-dimensional interpolation function for irregularly-spaced data.*
"""
@estimsolver InvDistWeight begin
  @param neighbors = nothing
  @param distance = Euclidean()
end

function solve(problem::EstimationProblem, solver::InvDistWeight)
  # retrieve problem info
  pdata = data(problem)
  pdomain = domain(problem)

  # result for each variable
  μs = []; σs = []

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # get variable type
      V = variables(problem)[var]

      # retrieve non-missing data
      locs = findall(!ismissing, pdata[var])
      𝒟 = view(pdata, locs)
      X = coordinates(𝒟)
      z = 𝒟[var]

      # number of data points for variable
      ndata = length(z)

      @assert ndata > 0 "estimation requires data"

      # allocate memory
      varμ = Vector{V}(undef, nelms(pdomain))
      varσ = Vector{V}(undef, nelms(pdomain))

      # fit search tree
      M = varparams.distance
      if M isa NearestNeighbors.MinkowskiMetric
        tree = KDTree(X, M)
      else
        tree = BruteTree(X, M)
      end

      # keep track of estimated locations
      estimated = falses(nelms(pdomain))

      # consider data locations as already estimated
      for (loc, datloc) in datamap(problem, var)
        estimated[loc] = true
        varμ[loc] = pdata[var][datloc]
        varσ[loc] = zero(V)
      end

      # determine number of nearest neighbors to use
      k = varparams.neighbors == nothing ? ndata : varparams.neighbors

      @assert k ≤ ndata "number of neighbors must be smaller or equal to number of data points"

      # pre-allocate memory for coordinates
      coords = MVector{ncoords(pdomain),coordtype(pdomain)}(undef)

      # estimation loop
      for location in traverse(pdomain, LinearPath())
        if !estimated[location]
          coordinates!(coords, pdomain, location)

          idxs, dists = knn(tree, coords, k)

          weights = one(V) ./ dists
          weights /= sum(weights)

          values = view(z, idxs)

          varμ[location] = sum(weights[i]*values[i] for i in eachindex(values))
          varσ[location] = minimum(dists)
        end
      end

      push!(μs, var => varμ)
      push!(σs, var => varσ)
    end
  end

  EstimationSolution(pdomain, Dict(μs), Dict(σs))
end
