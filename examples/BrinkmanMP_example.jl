using eFEMpart

function brinkmanMP_example()
  x0 = 2.0
  mesh = squareMeshFluid([-x0,x0,-1,1],32)

  μ0 = 1.0
  a10 = 0.0
  a20 = 1.0

  μ  = μ0*ones(length(mesh.xy))
  a1  = a10*ones(length(mesh.xy))
  a2  = a20*ones(length(mesh.xy)) 
  G  = 1e0
  param = BrinkmanMPParam(μ,a1,a2)

  OperatorType = :BrinkmanMP2D

  dNodes = Dirichlet(:left,:top,:bottom)
  nNodes = Neumann(:right)
  Nodes = [dNodes,nNodes]

  dUBCf = Dirichlet((x,y) -> (x==-x0 ? 1.0-y^2 : 0.0))
  dVBCf = Dirichlet((x,y) -> 0.0)
  bcfun = [dUBCf,dVBCf]

  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  sol = solve(prob,mesh,param)
  
  vtkname = Path("solutions","brinkmanMP_example")
  sd = ScalarData(sol.p)
  sn = ScalarNames("pressure")
  vd = VectorData([sol.u,sol.v])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  return nothing
end
 
