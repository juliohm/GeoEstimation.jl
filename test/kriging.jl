@testset "Kriging" begin
  # -----------------
  # 1D PROBLEM (ALL)
  # -----------------

  data1D = georef((z=[0.0,0.1,0.2,0.3,0.4,0.5,0.4,0.3,0.2,0.1,0.0],), collect(0.0:10.0:100.0)')
  grid1D = CartesianGrid(100)
  problem = EstimationProblem(data1D, grid1D, :z)

  global_kriging = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),))
  nearest_kriging = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.), maxneighbors=3))
  local_kriging = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),
                                 maxneighbors=3, neighborhood=NormBall(100.)))

  solvers   = [global_kriging, nearest_kriging, local_kriging]
  solnames  = ["global", "nearest", "local"]
  solutions = [solve(problem, solver) for solver in solvers]

  if visualtests
    for i in 1:3
      solution, sname = solutions[i], solnames[i]
      @test_reference "data/krig-1D-$(sname).png" plot(solution)
    end
  end

  # --------------------
  # 2D PROBLEM (GLOBAL)
  # --------------------

  data2D = georef((z=[1.,0.,1.],), [(25.,25.), (50.,75.), (75.,50.)])
  grid2D = CartesianGrid(100,100)

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),))

  solution = solve(problem, solver)

  # basic checks
  inds = LinearIndices(size(grid2D))
  S = solution[:z]
  @test isapprox(S[inds[26,26]], 1., atol=1e-3)
  @test isapprox(S[inds[51,76]], 0., atol=1e-3)
  @test isapprox(S[inds[76,51]], 1., atol=1e-3)

  if visualtests
    @test_reference "data/krig-2D-global.png" contourf(solution)
  end

  # ---------------------
  # 2D PROBLEM (NEAREST)
  # ---------------------

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.), maxneighbors=3))

  solution = solve(problem, solver)

  # basic checks
  inds = LinearIndices(size(grid2D))
  S = solution[:z]
  @test isapprox(S[inds[26,26]], 1., atol=1e-3)
  @test isapprox(S[inds[51,76]], 0., atol=1e-3)
  @test isapprox(S[inds[76,51]], 1., atol=1e-3)

  if visualtests
    @test_reference "data/krig-2D-nearest.png" contourf(solution)
  end

  # -------------------
  # 2D PROBLEM (LOCAL)
  # -------------------

  problem = EstimationProblem(data2D, grid2D, :z)

  solver = Kriging(:z => (variogram=GaussianVariogram(range=35.,nugget=0.),
                          maxneighbors=3, neighborhood=NormBall(100.)))

  solution = solve(problem, solver)

  # basic checks
  inds = LinearIndices(size(grid2D))
  S = solution[:z]
  @test isapprox(S[inds[26,26]], 1., atol=1e-3)
  @test isapprox(S[inds[51,76]], 0., atol=1e-3)
  @test isapprox(S[inds[76,51]], 1., atol=1e-3)

  if visualtests
    @test_reference "data/krig-2D-local.png" contourf(solution)
  end
end
