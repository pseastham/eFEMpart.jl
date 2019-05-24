using eFEMpart
include("Laplace_example_param.jl")

function laplace_example(p::Lap_Para)
  # define mesh through eFEMpart code
  mesh = squareMesh([p.xmin,p.xmax,p.ymin,p.ymax],p.Nintervals,p.basisorder)

  # define variable name and operator type
  OperatorType = :Poisson2D

  # Define auxiliary information
  dNodes = Dirichlet(:left,:right)
  nNodes = Neumann(:top,:bottom)

  # Define dirichlet, neumann, and forcing functions
  dBCf = Dirichlet((x,y) -> (x==p.xmin ? p.lBCval : p.rBCval))
  nBCf = Neumann((x,y) -> 0.0)
  ff   = Forcing((x,y) -> 0.0)

  # required technicality for grouping auxiliary information
  Nodes = [dNodes,nNodes]; bcfun = [dBCf,nBCf,ff]

  # define problem variable
  prob = Problem(mesh,Nodes,bcfun,OperatorType)

  # solve problem
  sol = solve(prob,mesh)

  # save solution
  vtkname = Path(p.solFolder,p.solFile)
  sd = ScalarData(sol.u)
  sn = ScalarNames(p.solName)
  vtksave(mesh,sd,sn,vtkname)
end

function main()
  p = Lap_Para()
  laplace_example(p)
end

main()
