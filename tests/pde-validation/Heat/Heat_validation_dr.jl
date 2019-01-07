# NOT CURRENTLY WORKING !!! NEED TO ADD TIME STEPPING ROUTINES TO eFEM

using eFEM, JLD

function main()
  function computeNorms(dt,i)
    tic()
    mesh = squareMesh([-2,2.0,-1.,1.],16,2);
    Tfinal = 0.5
    varname      = "temperature"
    OperatorType = :Poisson2D     # -Delta u = f

    # Define auxiliary information
    dNodes = Dirichlet(:left,:right)
    nNodes = Neumann(:bottom,:top)
    dBCf = Dirichlet((x,y) -> 0.0)
    nBCf = Neumann((x,y) -> 0.0)

    # required technicality for grouping auxiliary information
    Nodes = [dNodes,nNodes]; bcfun = [dBCf,nBCf]

    # define problem variable (helps solver)
    prob = Problem(mesh,Nodes,bcfun,varname,OperatorType)

    param = PoissonParam(1.0)
    xic = [mesh.xy[i].x for i=1:length(mesh.xy)]
    U0  = sin.(pi*xic/2)

    U = FEMForwardEuler(prob,mesh,param,U0,Tfinal,dt)
    sol = ScalarSolution(U)

    #fn = Path("HeatValidation$(i)")
    # print problem
    #vtksave(mesh,prob,param,sol,fn)

    time=toq()

    # generate exact solution
    Uexact(x,y) = exp(-(pi/2)^2*Tfinal)*sin(pi*x/2)
    UexactArr = [Uexact.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]

    # compute difference
    err = sol.u - UexactArr

    # compute Li norm
    L1   = DomainNorm(mesh.xy,mesh.cm,err;normID="1")
    L2   = DomainNorm(mesh.xy,mesh.cm,err;normID="2")
    Linf = DomainNorm(mesh.xy,mesh.cm,err;normID="Inf")

    return L1,L2,Linf,time
  end

  start = 1e-4

  dtarr    = [start,start,start/2,start/4,start/8,start/16]
  N       = length(dtarr)
  L1arr   = zeros(N)
  L2arr   = zeros(N)
  LInfarr = zeros(N)
  timearr = zeros(N)

  for i=1:N
      L1arr[i],L2arr[i],LInfarr[i],timearr[i] = computeNorms(dtarr[i],i)
      println("completed i=",i,"/",N)
  end

  L1arr = L1arr[2:end]
  L2arr = L2arr[2:end]
  LInfarr = LInfarr[2:end]
  timearr = timearr[2:end]
  dtarr = dtarr[2:end]


  run(`mkdir -p data`)

  save("data/TEMPHeat_validation.jld", "L1arr", L1arr,
                                        "L2arr", L2arr,
                                        "LInfarr",LInfarr,
                                        "timearr",timearr,
                                        "dtarr",dtarr)
  nothing
end

main()
