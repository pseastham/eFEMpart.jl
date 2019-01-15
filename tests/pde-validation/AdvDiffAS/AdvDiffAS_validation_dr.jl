using eFEMpart, JLD

using LinearAlgebra # for condition number calculation

function main()
  function computeNorms(N,order)
    start = time()
    # load mesh
    mesh = squareMesh([-2,2,0,1],N,order)
    
    varname      = "chemicalAxisym"
    OperatorType = :AdvDiffAS
   
    κ = 2.3
    windX(x,y) = -x^4*y^2 + 5.9; windY(x,y) = x^3*y^3
    wx = [windX.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    wy = [windY.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    param  = AdvDiffParam(wx,wy,κ)

    dNodes = Dirichlet(:left,:right,:top,:bottom)
    dBCf = Dirichlet((x,y) -> x^4*y^4)
    ff   = Forcing((x,y) -> 4*5.9*x^3*y^4 - 4*κ*x^2*y^2*(4*x^2 + 3*y^2))
    Nodes = [dNodes]
    bcfun = [dBCf,ff]

    prob = Problem(mesh,Nodes,bcfun,varname,OperatorType)
    sol = solve(prob,mesh,param)

    elapsed = time() - start
    
    # compute condition number of operator matrix
    #LinOp = GenerateSystem(mesh,prob,param)
    #ApplyBC!(LinOp,mesh,prob,param,OperatorType)
    #κ = cond(Array(LinOp.Op),2)
    κ = 0.5

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
    
    return κ,h,L1,L2,Linf,elapsed
  end

  Narr    = [8,8,16,32]#,64]#,128]#,256]
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

  save("data/TEMP_advdiffAS_validation.jld", "harr",harr[:,1],
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
