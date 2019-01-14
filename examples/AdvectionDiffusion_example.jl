using eFEMpart

function ad_example()
  Nintervals = 32
  basisorder=2
  xmin = -1.0; xmax = 1.0
  ymin = -1.0; ymax = 1.0
  square = [xmin,xmax,ymin,ymax]
  mesh = squareMesh(square,Nintervals,basisorder)

  varname      = "chemical"
  OperatorType = :AdvDiff2D
  
  dNodes = Dirichlet(:left,:right,:top,:bottom)
  dBCf   = Dirichlet((x,y) -> x^4*y^4)
  ff     = Forcing((x,y) -> 8*x^3*y^5*(4-x^2) - 8*x^5*y^3*(1-y^2) - 12*x^2*y^2*(x^2+y^2))
  
  # parameter
  κ = 1.0
  windX(x,y) =  2*y*(4-x^2);  windY(x,y) = -2*x*(1-y^2)
  wx = [windX.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
  wy = [windY.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
  param  = AdvDiffParam(wx,wy,κ)

  Nodes = [dNodes]
  bcfun = [dBCf,ff]
  
  prob = Problem(mesh,Nodes,bcfun,varname,OperatorType)
  sol = solve(prob,mesh,param)

  # save solution
  vtkname = Path("solutions","ad_example")
  sd = ScalarData(sol.u)
  sn = ScalarNames("chemical")
  vd = VectorData([wx,wy])
  vn = VectorNames("velocity")
  vtksave(mesh,sd,sn,vd,vn,vtkname)

  return nothing
end

