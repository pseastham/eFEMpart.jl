using eFEM

function main()
  tic()
  N = 128; order=2
  mesh = squareMesh([-2.,2.,-1.,1.],N,order)
  t1 = toq()

  tic()
  varname      = "temperature"
  OperatorType = :Poisson2D

  # Define auxiliary information
  dNodes = Dirichlet(:left,:right,:top)
  nNodes = Neumann(:bottom)
  dBCf = Dirichlet((x,y) -> 0.25*(2.0-x)^2)
  nBCf = Neumann((x,y) -> 0.0)
  ff   = Forcing((x,y) -> 0.0)

  # required technicality for grouping auxiliary information
  Nodes = [dNodes,nNodes]; bcfun = [dBCf,nBCf,ff]

  # define problem variable
  prob = Problem(mesh,Nodes,bcfun,varname,OperatorType)
  t2 = toq()

  tic()
  # solve problem
  sol = solve(prob,mesh)
  t3 = toq()

  tic()
  # save solution
  vtksave(mesh,prob,sol)
  t4 = toq()

  # print timing results
  r1=3; r2=5
  ttot = t1+t2+t3+t4
  t1frac = round(t1/ttot,r1)
  t2frac = round(t2/ttot,r1)
  t3frac = round(t3/ttot,r1)
  t4frac = round(t4/ttot,r1)
  println("N = ",N,"\t\t\t\t\t\t(%TOTAL)")
  println("mesh creation took .......... ",round(t1,r2)," seconds \t(",t1frac,"%)")
  println("problem definition took ..... ",round(t2,r2)," seconds \t(",t2frac,"%)")
  println("solution took ............... ",round(t3,r2)," seconds \t(",t3frac,"%)")
  println("vtk save took ............... ",round(t4,r2)," seconds \t(",t4frac,"%)")
end
