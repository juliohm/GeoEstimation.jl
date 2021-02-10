@testset "IDW" begin
  # basic problem
  geodata = georef((variable=[1.,0.,1.],), [25. 50. 75.;  25. 75. 50.])
  domain  = RegularGrid(100,100)
  problem = EstimationProblem(geodata, domain, :variable)

  solver = IDW(:variable => (neighbors=3,))

  solution = solve(problem, solver)

  if visualtests
    @test_ref_plot "data/idw.png" contourf(solution,size=(800,400))
  end

  # haversine distance
  geodata = georef((variable=[4.0,-1.0,3.0],), [50.0 100.0 200.0; -30.0 30.0 10.0])
  domain  = RegularGrid((1.0, -89.0), (359.0, 89.0), dims=(200, 100))
  problem = EstimationProblem(geodata, domain, :variable)

  solver = IDW(:variable => (neighbors=3, distance=Haversine(1.0)))

  solution = solve(problem, solver)

  if visualtests
    @test_ref_plot "data/idw-haversine.png" contourf(solution,size=(800,200))
  end
end
