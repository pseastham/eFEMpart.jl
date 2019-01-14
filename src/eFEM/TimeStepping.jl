# File to be loaded into eFEMpart

"""
FEMForwardEuler()

Executes Forward Euler time stepping to solve the linear equation Mx=b
and also computes the mass matrix M and b results from a matrix-vector product
using the linear operator. Uses LU factorization on the mass matrix
for the first step, then each step after that costs only O(n^2)

  Right now written to accept :Poisson2D OperatorType
"""
function FEMForwardEuler(prob,mesh,param,U0,Tfinal,dt;
                         printVTK=false,printTotal=100)
  t            = 0.0
  Ndiv = ceil(Tfinal/dt/printTotal)

  # create mass matrix
  M = MassMatrix(mesh)

  # perform LU decomp on mass matrix
  Mlu = lufact(M)

  # obtain weak form of initial condition
  U = forceassemble(mesh,U0)
  sol = ScalarSolution(U)

  # generate derivative matrices
  LinOp = MatrixGen(mesh,prob,param)

  # Apply boundary conditions
  BCApply!(LinOp,mesh,prob,param,prob.OperatorType)

  A = LinOp.Op
  F = LinOp.rhs

  i = 1; j=1
  # step in time
  while t<Tfinal
    # compute RHS
    RHS = (M-dt*A)*U + dt*F
    # use LU factorization to compute solution
    U = Mlu\RHS

    # print problem
    if printVTK
      if mod(j,Ndiv)==0
        s = string(prob.name)
        fn = Path(string(s,"_",i))
        sol.u = U
        vtksave(mesh,prob,param,sol,fn)
        i += 1
      end
    end
    # step time
    t += dt
    j += 1
  end

  return U
end

"""
  progressBar(t,timer1,timer2,timeAvg,counter,
              printTotal,printCounter,printSkip)

Calculates and displays progress bar information
"""
function progressBar(t,timer,counter,printTotal,printCounter,printSkip)
  timer[1] += toq(); tic()
  percentage = Int(floor(printCounter/printTotal*100))

  progressBarDisplay(percentage,printCounter,printTotal,timer)
end


"""
  progressBarDisplay()

Displays progress bar
"""
function progressBarDisplay(percentage,printCounter,printTotal,timer)
  const TERMINALLENGTH = 80

  percent = string("\r  [",percentage,"%]")
  skip1  = 14-length(percent)
  exportFigure = string("exported figure ",printCounter,"/",printTotal)
  skip2  = 42-(length(percent)+skip1+length(exportFigure))
  timeTaken = string("stopwatch: ",sToHMS(timer[1]),"\r")

  iterChar(TERMINALLENGTH," ")

  print(percent)

  space()
  iterChar(skip1,"-")
  space()
  print(exportFigure)

  space()
  iterChar(skip2,"-")
  space()
  print(timeTaken)

  if printCounter==printTotal
    println()
  end
end


"""
  iterChar(iterations::Int,character::String)

Prints string character iterations number of times. Used in progress bar
"""
function iterChar(iterations::Int,character::String)
  for i=1:iterations
    print(character)
  end
end

"""
  space()

Prints single space
"""
function space()
  print(" ")
end

"""
  scalarIC(xy,dNodes,dBCarr=0)

creates initial condition array for scalar equation using dirichlet Boundary nodes. If no dirichlet Boundary nodes are given, just returns zero vector
"""
function scalarIC(xy,dNodes,dBCarr)
  N = length(xy)
  ICarr = zeros(N)
  if length(dNodes) != 0
    ICarr[dNodes] = dBCarr
  end

  return ICarr
end

scalarIC(xy) = zeros(Float64,length(xy))

"""
  sToHMS(s::Float64)

Converts Float64 in seconds to HH:MM:SS string format
"""
function sToHMS(s)
  minutes = 0.0
  hours   = 0.0
  timestr = ""

  # compute times
  hours   = Int(floor(s/3600))
  minutes = Int(floor(s/60 - 60*hours))
  seconds = Int(round(s - 3600*hours - 60*minutes,0))

  hours < 10   ? hstr = "0" : hstr = ""
  minutes < 10 ? mstr = "0" : mstr = ""
  seconds < 10 ? sstr = "0" : sstr = ""

  if hours > 0
    timestr=string(timestr,hstr,hours,"h:")
  end
  if (minutes>0) || (hours > 0)
    timestr=string(timestr,mstr,minutes,"m:")
  end

  timestr = string(timestr,sstr,seconds,"s")

  return timestr
end
