using eFEMpart
include("Stokes_example_param.jl")

function stokes_example(p::Stokes_Para)
  mesh = squareMeshFluid([p.xmin,p.xmax,p.ymin,p.ymax],p.Nintervals)

  param = FluidParam(p.μ)
  OperatorType = :Stokes2D

  dNodes = Dirichlet(:left,:bottom,:top)
  nNodes = Neumann(:right)
  Nodes = [dNodes,nNodes]

  dUBCf = Dirichlet((x,y) -> (x==p.xmin ? -p.lBCval*(p.ymax-y)*(p.ymin-y) : 0.0))
  dVBCf = Dirichlet((x,y) -> 0.0)
  bcfun = [dUBCf,dVBCf]

  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  sol = solve(prob,mesh,param)

  # save solution
  vtkname = Path(p.solFolder,p.solFile)
  sd = ScalarData(sol.p)
  sn = ScalarNames("pressure")
  vd = VectorData([sol.u,sol.v])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  nothing
end

function main()
  p = Stokes_Para()
  stokes_example(p)
end

main()
