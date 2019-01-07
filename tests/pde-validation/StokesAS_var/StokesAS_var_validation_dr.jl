using eFEM, JLD

function main()
  function computeNorms(N)
    tic()
    # load mesh
    mesh = squareMeshFluid([-2,2,0,1],N)

    μ0     = 2.3
    μ(x,y) = μ0 + exp(x*y)
    xm = [i.x for i in mesh.xy]
    ym = [i.y for i in mesh.xy]
    μarr = [μ(xm[i],ym[i]) for i=1:length(mesh.xy)]
    param  = FluidParam(μarr)

    OperatorType = :StokesAS

    dUNodes = Dirichlet(:left,:top,:right)
    dVNodes = Dirichlet(:left,:top,:bottom,:right)
    Nodes = [dUNodes,dVNodes]

    dUBCf = Dirichlet((x,y) -> -x^4*y^2 + 5.7)
    dVBCf = Dirichlet((x,y) ->  x^3*y^3)

    dudx(x,y) = -4*x^3*y^2; dudy(x,y) = -2*x^4*y
    dvdx(x,y) = 3*x^2*y^3;  dvdy(x,y) = 3*x^3*y^2
    cross(x,y) = dudy(x,y) + dvdx(x,y)

    d2dudx(x,y) = -12*x^2*y^2; d2dudy(x,y) = -2*x^4
    d2dvdx(x,y) = 6*x*y^3;     d2dvdy(x,y) = 6*x^3*y
    #lapu(x,y) = d2dudx(x,y) + d2dudy(x,y) + dudy(x,y)/y
    lapu(x,y) = d2dudx(x,y) + d2dudy(x,y) - 2*x^4
    #lapvplus(x,y) = d2dvdx(x,y) + d2dvdy(x,y) + dvdy(x,y)/y - v/y^2  
    lapvplus(x,y) = d2dvdx(x,y) + d2dvdy(x,y) + 3*x^3*y  - x^3*y

    dpdx(x,y) = 3*x^2*y^3; dpdy(x,y) = 3*x^3*y^2

    dμdx(x,y) = y*exp(x*y)
    dμdy(x,y) = x*exp(x*y)
    
    Fx(x,y) = -μ(x,y)*lapu(x,y) - 2*dμdx(x,y)*dudx(x,y) - dμdy(x,y)*cross(x,y) + dpdx(x,y)
    Fy(x,y) = -μ(x,y)*lapvplus(x,y) - 2*dμdy(x,y)*dvdy(x,y) - dμdx(x,y)*cross(x,y) + dpdy(x,y)

    ffx   = Forcing(Fx)
    ffy   = Forcing(Fy)
    bcfun = [dUBCf,dVBCf,ffx,ffy]

    prob = Problem(mesh,Nodes,bcfun,OperatorType)
    sol = solve(prob,mesh,param)

    time=toq()

    # compute condition number of operator matrix
    #LinOp = GenerateSystem(mesh,prob,param)
    #ApplyBC!(LinOp,mesh,prob,param,OperatorType)
    #κ = cond(Array(LinOp.Op),2)
    κ = 0.5

    # compute mesh h
    h = hCalc(mesh)

    # generate exact solution
    Uexact(x,y) = dUBCf.f(x,y)
    Vexact(x,y) = dVBCf.f(x,y)
    Pexact(x,y) = x^3*y^3
    UexactArr = [Uexact.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    VexactArr = [Vexact.(mesh.xy[i].x,mesh.xy[i].y) for i=1:length(mesh.xy)]
    PexactArr = [Pexact.(mesh.xyp[i].x,mesh.xyp[i].y) for i=1:length(mesh.xyp)]

    # plot solution
    sD = ScalarData(sol.p)
    sN = ScalarNames("pressure")
    vD = VectorData([sol.u,sol.v])
    vN = VectorNames("velocity")
    fn = Path("solution_var")
    vtksave(mesh,sD,sN,vD,vN,fn)

    # compute difference
    velErr = sqrt.((sol.u - UexactArr).^2 + (sol.v -VexactArr).^2)
    presErr = sol.p - PexactArr

    # compute Li norm
    L1v   = DomainNorm(mesh.xy,mesh.cm,velErr;normID="1")
    L2v   = DomainNorm(mesh.xy,mesh.cm,velErr;normID="2")
    Linfv = DomainNorm(mesh.xy,mesh.cm,velErr;normID="Inf")
    L1p   = DomainNorm(mesh.xyp,mesh.cmp,presErr;normID="1")
    L2p   = DomainNorm(mesh.xyp,mesh.cmp,presErr;normID="2")
    Linfp = DomainNorm(mesh.xyp,mesh.cmp,presErr;normID="Inf")

    return κ,h,L1v,L2v,Linfv,L1p,L2p,Linfp,time
  end

  Narr     = [4,4,8,16,32]#,64]#,96]
  N        = length(Narr)
  harr     = zeros(N)
  L1arrV   = zeros(N)
  L2arrV   = zeros(N)
  LInfarrV = zeros(N)
  L1arrP   = zeros(N)
  L2arrP   = zeros(N)
  LInfarrP = zeros(N)
  timearr  = zeros(N)
  κarr     = zeros(N)

  for i=1:N
    n = Narr[i]
    κarr[i],harr[i],L1arrV[i],L2arrV[i],LInfarrV[i],
      L1arrP[i],L2arrP[i],LInfarrP[i],timearr[i] =
      computeNorms(n)
      println("completed N=$(n)")
  end

  κarr     = κarr[2:end]
  harr     = harr[2:end]
  L1arrV   = L1arrV[2:end]
  L2arrV   = L2arrV[2:end]
  LInfarrV = LInfarrV[2:end]
  L1arrP   = L1arrP[2:end]
  L2arrP   = L2arrP[2:end]
  LInfarrP = LInfarrP[2:end]
  timearr  = timearr[2:end]

  # export output to *.jld file
  run(`mkdir -p data`)

  save("data/TEMP_stokesAS_var_validation.jld", "harr",harr,
                                          "VL1arr", L1arrV,
                                          "VL2arr", L2arrV,
                                          "VLInfarr",LInfarrV,
                                          "PL1arr", L1arrP,
                                          "PL2arr", L2arrP,
                                          "PLInfarr",LInfarrP,
                                          "timearr",timearr,
                                          "κarr",κarr)

  nothing
end

main()
