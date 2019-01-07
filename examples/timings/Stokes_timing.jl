using eFEM

# :Poisson2D        -> -Lap(u) = f
# :AdvDiff2D        -> w*Grad(u) - Lap(u) = f
# :Stokes2D         -> -Lap(u) + Grad(p) = F, Div(u) = 0
# :AdvDiffAxisym    -> AdvDiff in axisymmetric coordinates
# :StokesAxisym     -> Stokes in axisymmetric coordinates

# I'm pretty sure there's something wrong with the way forcing
#   functions are evaluated!! Test this! (WeakScalar in Assembly.jl)
function main()
  tic()
  N = 128
  mesh = squareMeshFluid([-2,2,-1,1],N)
  t1 = toq()

  tic()
  μ = 1.0
  param = FluidParam(μ)

  OperatorType = :Stokes2D

  dNodes = Dirichlet(:left,:top,:bottom)
  nNodes = Neumann(:right)
  Nodes = [dNodes,nNodes]

  dUBCf = Dirichlet((x,y) -> (x==-2.0 ? 1.0-y^2 : 0.0))
  dVBCf = Dirichlet((x,y) -> 0.0)
  bcfun = [dUBCf,dVBCf]

  prob = Problem(mesh,Nodes,bcfun,OperatorType)
  t2 = toq()

  tic()
  sol = solve(prob,mesh,param)
  t3 = toq()

  tic()
  vtksave(mesh,prob,param,sol)
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
