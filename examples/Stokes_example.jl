using eFEMpart

function stokes_example()
  mesh = squareMeshFluid([-2,2,-1,1],32)

  μ = 1.0
  param = FluidParam(μ)
  OperatorType = :Stokes2D

  dNodes = Dirichlet(:left,:bottom,:top)
  nNodes = Neumann(:right)
  Nodes = [dNodes,nNodes]

  dUBCf = Dirichlet((x,y) -> (x==-2.0 ? 1.0-y^2 : 0.0))
  dVBCf = Dirichlet((x,y) -> 0.0)
  bcfun = [dUBCf,dVBCf]

  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  sol = solve(prob,mesh,param)

  vtkname = Path("solutions","stokes_example")
  sd = ScalarData(sol.p)
  sn = ScalarNames("pressure")
  vd = VectorData([sol.u,sol.v])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

end
