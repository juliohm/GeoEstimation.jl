# GeoEstimation.jl

[![][build-img]][build-url] [![][codecov-img]][codecov-url]

Geostatistical estimation solvers for the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

### IDW

This solver provides a high-performance implementation of the inverse distance weighting scheme
introduced in the very early days of geostatistics (see [Shepard 1968](https://dl.acm.org/citation.cfm?id=810616)).
It is perhaps the simplest first attempt in the literature to perform estimation based on the
notion of proximity to data locations.

This implementation makes use of k-d trees from the [NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl)
package, which leads to a fast estimation method for large or high-resolution spatial domains.
Although this method is recommended for fast assessment of a new field, it has poor statistical
properties (lacks covariance model) and should mainly be used for qualitative purposes.

### LWR

This solver provides an implementation of locally weighted regression (a.k.a. LOESS) introduced by
[Cleveland 1979](http://www.stat.washington.edu/courses/stat527/s13/readings/Cleveland_JASA_1979.pdf).
It is the most natural generalization of inverse distance weighting in which one is allowed to use a
custom weight function instead of distance-based weights.

Like in the inverse distance weighting solver, this solver makes use of k-d trees from the
[NearestNeighbors.jl](https://github.com/KristofferC/NearestNeighbors.jl) package for fast data
lookup. Locally weighted regression (LWR or LOESS) is a popular non-parametric method, however
it still has poor statistical properties compared to other estimation methods such as Kriging
that explicitly model spatial correlation.

In the current implementation, the estimation variance is computed assuming Gaussian residuals.

### Kriging

This polyalgorithm solver provides an implementation of various forms of Kriging introduced by
[Matheron 1971](https://books.google.com.br/books/about/The_Theory_of_Regionalized_Variables_and.html).
Kriging is a popular method in various industries due to its good statistical properties and flexibility.
Unlike the previous solvers, Kriging relies on the specification of a variogram model, which can be
fit to geospatial data.

## Installation

Get the latest stable release with Julia's package manager:

```julia
] add GeoEstimation
```

## Usage

This package is part of the [GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework.

For a simple example of usage, please check the main documentation.

## Asking for help

If you have any questions, please contact our community on the [gitter channel](https://gitter.im/JuliaEarth/GeoStats.jl).

[build-img]: https://img.shields.io/github/workflow/status/JuliaEarth/GeoEstimation.jl/CI?style=flat-square
[build-url]: https://github.com/JuliaEarth/GeoEstimation.jl/actions

[codecov-img]: https://codecov.io/gh/JuliaEarth/GeoEstimation.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaEarth/GeoEstimation.jl
