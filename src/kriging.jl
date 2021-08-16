# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Kriging(var₁=>param₁, var₂=>param₂, ...)

A polyalgorithm Kriging estimation solver.

Each pair `var=>param` specifies the `KrigingParam` `param`
for the Kriging variable `var`. In order to avoid boilerplate
code, the constructor expects pairs of `Symbol` and `NamedTuple`
instead.

## Parameters

* `variogram` - Variogram model (default to `GaussianVariogram()`)
* `mean`      - Simple Kriging mean
* `degree`    - Universal Kriging degree
* `drifts`    - External Drift Kriging drift functions

Latter options override former options. For example, by specifying
`drifts`, the user is telling the algorithm to ignore `degree` and
`mean`. If no option is specified, Ordinary Kriging is used by
default with the `variogram` only.

* `minneighbors` - Minimum number of neighbors (default to `1`)
* `maxneighbors` - Maximum number of neighbors (default to `10`)
* `neighborhood` - Search neighborhood (default to `nothing`)
* `distance`     - Distance used to find nearest neighbors (default to `Euclidean()`)

The `maxneighbors` option can be used to perform approximate Kriging
with a subset of data points per estimation location. If `maxneighbors`
is set to `nothing`, global Kriging is performed instead. Two neighborhood
search methods are available depending on the value of `neighborhood`:

* If a `neighborhood` is provided, local Kriging is performed by sliding
  the `neighborhood` in the domain.

* If `neighborhood` is not provided, the Kriging system is built using
  `maxneighbors` nearest neighbors according to a `distance`.

## Examples

Solve the variable `:var₁` with Simple Kriging by specifying
the `mean`, and the variable `:var₂` with Universal Kriging
by specifying the `degree` and the `variogram` model.

```julia
julia> Kriging(
  :var₁ => (mean=1.,),
  :var₂ => (degree=1, variogram=SphericalVariogram(range=20.))
)
```

Solve all variables of the problem with the default parameters
(i.e. Ordinary Kriging with unit Gaussian variogram):

```julia
julia> Kriging()
```
"""
@estimsolver Kriging begin
  @param variogram = GaussianVariogram()
  @param mean = nothing
  @param degree = nothing
  @param drifts = nothing
  @param minneighbors = 1
  @param maxneighbors = 10
  @param neighborhood = nothing
  @param distance = Euclidean()
end

function preprocess(problem::EstimationProblem, solver::Kriging)
  # retrieve problem info
  pdomain = domain(problem)
  pdata   = data(problem)
  ndata   = nelements(pdata)

  # result of preprocessing
  preproc = Dict{Symbol,NamedTuple}()

  for covars in covariables(problem, solver)
    for var in covars.names
      # get user parameters
      varparams = covars.params[(var,)]

      # determine which Kriging variant to use
      if varparams.drifts ≠ nothing
        estimator = ExternalDriftKriging(varparams.variogram, varparams.drifts)
      elseif varparams.degree ≠ nothing
        estimator = UniversalKriging(varparams.variogram, varparams.degree, embeddim(pdomain))
      elseif varparams.mean ≠ nothing
        estimator = SimpleKriging(varparams.variogram, varparams.mean)
      else
        estimator = OrdinaryKriging(varparams.variogram)
      end

      # determine minimum/maximum number of neighbors
      minneighbors = varparams.minneighbors
      maxneighbors = varparams.maxneighbors

      # determine neighborhood search method
      if varparams.maxneighbors ≠ nothing
        # upper bound in maxneighbors
        maxneighbors > ndata && (maxneighbors = ndata)
        if varparams.neighborhood ≠ nothing
          # local search with a neighborhood
          neigh = varparams.neighborhood

          if neigh isa MetricBall
            bsearcher = KBallSearch(pdata, maxneighbors, neigh)
          else
            searcher  = BallSearch(pdata, neigh)
            bsearcher = BoundedSearch(searcher, maxneighbors)
          end
        else
          # nearest neighbor search with a distance
          distance  = varparams.distance
          bsearcher = KNearestSearch(pdata, maxneighbors, metric=distance)
        end
      else
        # use all data points as neighbors
        bsearcher = nothing
      end

      # save preprocessed input
      preproc[var] = (estimator=estimator,
                      minneighbors=minneighbors,
                      maxneighbors=maxneighbors,
                      bsearcher=bsearcher)
    end
  end

  preproc
end

function solve(problem::EstimationProblem, solver::Kriging)
  # preprocess user input
  preproc = preprocess(problem, solver)

  # results for each variable
  μs = []; σs = []
  for var in name.(variables(problem))
    if preproc[var].maxneighbors ≠ nothing
      # perform Kriging with reduced number of neighbors
      varμ, varσ = solve_approx(problem, var, preproc)
    else
      # perform Kriging with all data points as neighbors
      varμ, varσ = solve_exact(problem, var, preproc)
    end

    push!(μs, var => varμ)
    push!(σs, Symbol(var,"_variance") => varσ)
  end

  georef((; μs..., σs...), domain(problem))
end

function solve_approx(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)

    mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

    # unpack preprocessed parameters
    estimator, minneighbors, maxneighbors, bsearcher = preproc[var]

    # determine value type
    V = mactypeof[var]

    # pre-allocate memory for result
    varμ = Vector{V}(undef, nelements(pdomain))
    varσ = Vector{V}(undef, nelements(pdomain))

    # pre-allocate memory for neighbors
    neighbors = Vector{Int}(undef, maxneighbors)

    # estimation loop
    for location in traverse(pdomain, LinearPath())
      pₒ = centroid(pdomain, location)

      # find neighbors with previously estimated values
      nneigh = search!(neighbors, pₒ, bsearcher)

      # skip location in there are too few neighbors
      if nneigh < minneighbors
        varμ[location] = NaN
        varσ[location] = NaN
      else
        # final set of neighbors
        nview = view(neighbors, 1:nneigh)

        # view neighborhood with data
        𝒟 = view(pdata, nview)

        # fit estimator to data
        krig = fit(estimator, 𝒟, var)

        # retrieve element at location
        uₒ = pdomain[location]

        # save mean and variance
        μ, σ² = predict(krig, uₒ)

        varμ[location] = μ
        varσ[location] = σ²
      end
    end

    varμ, varσ
end

function solve_exact(problem::EstimationProblem, var::Symbol, preproc)
    # retrieve problem info
    pdata = data(problem)
    pdomain = domain(problem)

    mactypeof = Dict(name(v) => mactype(v) for v in variables(problem))

    # unpack preprocessed parameters
    estimator, minneighbors, maxneighbors, bsearcher = preproc[var]

    # determine value type
    V = mactypeof[var]

    # pre-allocate memory for result
    varμ = Vector{V}(undef, nelements(pdomain))
    varσ = Vector{V}(undef, nelements(pdomain))

    # retrieve non-missing data
    locs = findall(!ismissing, pdata[var])
    𝒟 = view(pdata, locs)

    # fit estimator once
    krig = fit(estimator, 𝒟, var)

    # predict at all locations
    for location in traverse(pdomain, LinearPath())
      uₒ = pdomain[location]

      μ, σ² = predict(krig, uₒ)

      varμ[location] = μ
      varσ[location] = σ²
    end

    varμ, varσ
end
