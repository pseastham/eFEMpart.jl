using eFEMpart

function laplace_example()
  # define mesh through eFEMpart code; 
  # alternative is to load GMSH mesh
  Nintervals = 32
  basisorder = 2
  xmin = -2.0; xmax = 2.0
  ymin = -1.0; ymax = 1.0
  mesh = squareMesh([xmin,xmax,ymin,ymax],Nintervals,basisorder)

  # define variable name and operator type
  varname      = "temperature"
  OperatorType = :Poisson2D

  # Define auxiliary information
  dNodes = Dirichlet(:left,:right)
  nNodes = Neumann(:top,:bottom)

  # Define dirichlet, neumann, and forcing functions
  dBCf = Dirichlet((x,y) -> (x==-2.0 ? 1.0 : 0.0))
  nBCf = Neumann((x,y) -> 0.0)
  ff   = Forcing((x,y) -> 0.0)

  # required technicality for grouping auxiliary information
  Nodes = [dNodes,nNodes]; bcfun = [dBCf,nBCf,ff]

  # define problem variable
  prob = Problem(mesh,Nodes,bcfun,varname,OperatorType)

  # solve problem
  sol = solve(prob,mesh)

  # save solution
  vtksave(mesh,prob,sol)
end
