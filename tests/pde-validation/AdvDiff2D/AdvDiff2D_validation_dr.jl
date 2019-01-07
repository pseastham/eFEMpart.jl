using eFEM, JLD

function main()
  function computeNorms(N,order)
    tic()
    # load mesh
    mesh = squareMesh([-2,2,-1,1],N,order)
    
    # define parameter  
    windX(x,y) =  2*y*(4-x^2);    windY(x,y) = -2*x*(1-y^2)
    wx = [windX.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    wy = [windY.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    κ  = 2.3
    param  = AdvDiffParam(wx,wy,κ)
    varname = "chemical"
    OperatorType = :AdvDiff2D
    
    dNodes = Dirichlet(:left,:top,:bottom,:right)  
    dBCf = Dirichlet((x,y) -> x^4*y^4)
    ff   = Forcing((x,y) -> 8*x^3*y^5*(4-x^2) - 8*x^5*y^3*(1-y^2) - κ*12*x^2*y^2*(x^2+y^2))
    
    Nodes = [dNodes]; bcfun = [dBCf,ff]
    prob = Problem(mesh,Nodes,bcfun,varname,OperatorType)

    sol = solve(prob,mesh,param)
    
    time=toq() 

    # compute condition number of operator matrix
    LinOp = GenerateSystem(mesh,prob,param)
    ApplyBC!(LinOp,mesh,prob,param,OperatorType)
    κ = cond(Array(LinOp.Op),2)
    #κ = -1.0

    # compute mesh h
    h = hCalc(mesh)
    
    # generate exact solution
    Uexact(x,y) = dBCf.f(x,y)
    UexactArr = [Uexact.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    
    # compute difference
    err = sol.u - UexactArr
    
    # compute Li norm
    L1   = DomainNorm(mesh.xy,mesh.cm,err;normID="1")
    L2   = DomainNorm(mesh.xy,mesh.cm,err;normID="2")
    Linf = DomainNorm(mesh.xy,mesh.cm,err;normID="Inf")  
    
    return κ,h,L1,L2,Linf,time
  end

  Narr    = [4,4,8,16,32]#,64,128]#,256]
  order   = [1,2]
  N       = length(Narr)
  κarr    = zeros(N,2)
  harr    = zeros(N,2)
  L1arr   = zeros(N,2)
  L2arr   = zeros(N,2)
  LInfarr = zeros(N,2)
  timearr = zeros(N,2)

  for i=1:N
    for j=1:2
      n = Narr[i]
      o = order[j]
      κarr[i,j],harr[i,j],L1arr[i,j],L2arr[i,j],LInfarr[i,j],timearr[i,j] = computeNorms(n,o)
      println("completed N=$(n), order=$(o)")
    end
  end

  κarr = κarr[2:end,:]
  harr = harr[2:end,:]
  L1arr = L1arr[2:end,:]
  L2arr = L2arr[2:end,:]
  LInfarr = LInfarr[2:end,:]
  timearr = timearr[2:end,:]

  # export output to *.jld file
  run(`mkdir -p data`)

  save("data/TEMPadvdiff2D_validation.jld", "harr",harr[:,1],
                                          "o1L1arr", L1arr[:,1],
                                          "o1L2arr", L2arr[:,1],
                                          "o1LInfarr",LInfarr[:,1],
                                          "o1timearr",timearr[:,1],
                                          "o2L1arr", L1arr[:,2],
                                          "o2L2arr", L2arr[:,2],
                                          "o2LInfarr",LInfarr[:,2],
                                          "o2timearr",timearr[:,2],
                                          "κarr",κarr[:,2])

  nothing
end

main()
