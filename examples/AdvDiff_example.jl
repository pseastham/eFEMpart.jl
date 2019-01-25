using eFEMpart
include("AdvDiff_example_param.jl")

function ad_example(p::AD_Para)
  square = [p.xmin,p.xmax,p.ymin,p.ymax]
  mesh = squareMesh(square,p.Nintervals,p.basisorder)

  varname      = ""
  OperatorType = :AdvDiff2D
  
  dNodes = Dirichlet(:left,:top)
  nNodes = Neumann(:right,:bottom)

  dBCf   = Dirichlet((x,y) -> (x==p.xmin ? p.dirVal : 0.0))
  ff     = Forcing((x,y) -> 0.0)
  
  # parameter
  κ = 1/p.Pe
  windX(x,y) =  -y;  windY(x,y) = x
  wx = [windX.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
  wy = [windY.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
  param  = AdvDiffParam(wx,wy,κ)

  Nodes = [dNodes]
  bcfun = [dBCf,ff]
  
  prob = Problem(mesh,Nodes,bcfun,varname,OperatorType)
  sol = solve(prob,mesh,param)

  # save solution
  vtkname = Path(p.solFolder,p.solFile)
  sd = ScalarData(sol.u)
  sn = ScalarNames(p.solName)
  vd = VectorData([wx,wy])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  return nothing
end

function main()
  p = AD_Para()
  ad_example(p)
end

main()
