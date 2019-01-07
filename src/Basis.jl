# File loaded into eFEM

"""
  ~~~  EFFICIENT ~~~

  shape2D(s::Float64,t::Float64,order::Int)

Computes shape functions of order order (e.g. 1=linear, 2=quadratic)

ASSUMPTIONS:
- Only available for quadratic meshes, and linear order (order=1)

INPUT:
  s:      [Float64]
             horizontal location in computation domain [-1,1]^2
  t:      [Float64]
             vertical location in computational domain [-1,1]^2
  order:  [Int]
             Order of basis function (linear or quadratic)

OUTPUT:
  phi:    [Vector{Float64}]
             Nodal basis functions evaluated at (s,t)
  dphids  [Vector{Float64}]
             Derivative in s-direction of nodal basis function
               evaluated at (s,t)
  dphidt:  [Vector{Float64}]
             Derivative in t-direction of nodal basis function
               evaluated at (s,t)
"""
function shape2D(s::Float64,t::Float64,order::Int)
  if order == 1
    phi    = Array{Float64}(4)
    dphids = Array{Float64}(4)
    dphidt = Array{Float64}(4)

    phi[1]    = 0.25*(1-s)*(1-t)
    phi[2]    = 0.25*(1+s)*(1-t)
    phi[3]    = 0.25*(1+s)*(1+t)
    phi[4]    = 0.25*(1-s)*(1+t)
    dphids[1] = -0.25*(1-t)
    dphids[2] =  0.25*(1-t)
    dphids[3] =  0.25*(1+t)
    dphids[4] = -0.25*(1+t)
    dphidt[1] = -0.25*(1-s)
    dphidt[2] = -0.25*(1+s)
    dphidt[3] =  0.25*(1+s)
    dphidt[4] =  0.25*(1-s)
  elseif order == 2
    phi    = Array{Float64}(9)
    dphids = Array{Float64}(9)
    dphidt = Array{Float64}(9)

    ells1  = 0.5*s*(s-1); ellt1 = 0.5*t*(t-1);
    ells2  = 1-(s*s);     ellt2 = 1-(t*t);
    ells3  = 0.5*s*(s+1); ellt3 = 0.5*t*(t+1);
    dells1 = s-0.5;       dellt1 = t-0.5;
    dells2 = -2*s;        dellt2 = -2*t;
    dells3 = s+0.5;       dellt3 = t+0.5;

    phi[1] = ells1*ellt1
    phi[2] = ells3*ellt1
    phi[3] = ells3*ellt3
    phi[4] = ells1*ellt3
    phi[5] = ells2*ellt1
    phi[6] = ells3*ellt2
    phi[7] = ells2*ellt3
    phi[8] = ells1*ellt2
    phi[9] = ells2*ellt2
    dphids[1] = dells1*ellt1
    dphids[2] = dells3*ellt1
    dphids[3] = dells3*ellt3
    dphids[4] = dells1*ellt3
    dphids[5] = dells2*ellt1
    dphids[6] = dells3*ellt2
    dphids[7] = dells2*ellt3
    dphids[8] = dells1*ellt2
    dphids[9] = dells2*ellt2
    dphidt[1] = ells1*dellt1
    dphidt[2] = ells3*dellt1
    dphidt[3] = ells3*dellt3
    dphidt[4] = ells1*dellt3
    dphidt[5] = ells2*dellt1
    dphidt[6] = ells3*dellt2
    dphidt[7] = ells2*dellt3
    dphidt[8] = ells1*dellt2
    dphidt[9] = ells2*dellt2
  end

  return phi,dphids,dphidt
end

"""
  ~~~  EFFICIENT ~~~
"""
function shape2D!(phi,dphids,dphidt,
                 s::Float64,t::Float64,order::Int)
  if order == 1
    phi[1]    = 0.25*(1-s)*(1-t)
    phi[2]    = 0.25*(1+s)*(1-t)
    phi[3]    = 0.25*(1+s)*(1+t)
    phi[4]    = 0.25*(1-s)*(1+t)
    dphids[1] = -0.25*(1-t)
    dphids[2] =  0.25*(1-t)
    dphids[3] =  0.25*(1+t)
    dphids[4] = -0.25*(1+t)
    dphidt[1] = -0.25*(1-s)
    dphidt[2] = -0.25*(1+s)
    dphidt[3] =  0.25*(1+s)
    dphidt[4] =  0.25*(1-s)
  elseif order == 2

    ells1  = 0.5*s*(s-1);  ellt1 = 0.5*t*(t-1);
    ells2  = 1-(s*s);      ellt2 = 1-(t*t);
    ells3  = 0.5*s*(s+1);  ellt3 = 0.5*t*(t+1);
    dells1 = s-0.5;       dellt1 = t-0.5;
    dells2 = -2*s;        dellt2 = -2*t;
    dells3 = s+0.5;       dellt3 = t+0.5;

    phi[1] = ells1*ellt1
    phi[2] = ells3*ellt1
    phi[3] = ells3*ellt3
    phi[4] = ells1*ellt3
    phi[5] = ells2*ellt1
    phi[6] = ells3*ellt2
    phi[7] = ells2*ellt3
    phi[8] = ells1*ellt2
    phi[9] = ells2*ellt2
    dphids[1] = dells1*ellt1
    dphids[2] = dells3*ellt1
    dphids[3] = dells3*ellt3
    dphids[4] = dells1*ellt3
    dphids[5] = dells2*ellt1
    dphids[6] = dells3*ellt2
    dphids[7] = dells2*ellt3
    dphids[8] = dells1*ellt2
    dphids[9] = dells2*ellt2
    dphidt[1] = ells1*dellt1
    dphidt[2] = ells3*dellt1
    dphidt[3] = ells3*dellt3
    dphidt[4] = ells1*dellt3
    dphidt[5] = ells2*dellt1
    dphidt[6] = ells3*dellt2
    dphidt[7] = ells2*dellt3
    dphidt[8] = ells1*dellt2
    dphidt[9] = ells2*dellt2
  end

  nothing
end

"""
  ~~~  EFFICIENT ~~~
"""
function shape2D!(phi,dphids,dphidt,s::Float64,t::Float64,order::Int,gpt::Int)
  if order == 1
    phi[1,gpt]    = 0.25*(1-s)*(1-t)
    phi[2,gpt]    = 0.25*(1+s)*(1-t)
    phi[3,gpt]    = 0.25*(1+s)*(1+t)
    phi[4,gpt]    = 0.25*(1-s)*(1+t)
    dphids[1,gpt] = -0.25*(1-t)
    dphids[2,gpt] =  0.25*(1-t)
    dphids[3,gpt] =  0.25*(1+t)
    dphids[4,gpt] = -0.25*(1+t)
    dphidt[1,gpt] = -0.25*(1-s)
    dphidt[2,gpt] = -0.25*(1+s)
    dphidt[3,gpt] =  0.25*(1+s)
    dphidt[4,gpt] =  0.25*(1-s)
  elseif order == 2
    ells1  = 0.5*s*(s-1);  ellt1 = 0.5*t*(t-1);
    ells2  = 1-(s*s);      ellt2 = 1-(t*t);
    ells3  = 0.5*s*(s+1);  ellt3 = 0.5*t*(t+1);
    dells1 = s-0.5;       dellt1 = t-0.5;
    dells2 = -2*s;        dellt2 = -2*t;
    dells3 = s+0.5;       dellt3 = t+0.5;

    phi[1,gpt] = ells1*ellt1
    phi[2,gpt] = ells3*ellt1
    phi[3,gpt] = ells3*ellt3
    phi[4,gpt] = ells1*ellt3
    phi[5,gpt] = ells2*ellt1
    phi[6,gpt] = ells3*ellt2
    phi[7,gpt] = ells2*ellt3
    phi[8,gpt] = ells1*ellt2
    phi[9,gpt] = ells2*ellt2
    dphids[1,gpt] = dells1*ellt1
    dphids[2,gpt] = dells3*ellt1
    dphids[3,gpt] = dells3*ellt3
    dphids[4,gpt] = dells1*ellt3
    dphids[5,gpt] = dells2*ellt1
    dphids[6,gpt] = dells3*ellt2
    dphids[7,gpt] = dells2*ellt3
    dphids[8,gpt] = dells1*ellt2
    dphids[9,gpt] = dells2*ellt2
    dphidt[1,gpt] = ells1*dellt1
    dphidt[2,gpt] = ells3*dellt1
    dphidt[3,gpt] = ells3*dellt3
    dphidt[4,gpt] = ells1*dellt3
    dphidt[5,gpt] = ells2*dellt1
    dphidt[6,gpt] = ells3*dellt2
    dphidt[7,gpt] = ells2*dellt3
    dphidt[8,gpt] = ells1*dellt2
    dphidt[9,gpt] = ells2*dellt2
  end

  nothing
end

"""
  ~~~  EFFICIENT ~~~
  shapeEval(coef::Vector{FLoat64},basis::Vector{Float64})

Evaluates the sum given by the basis values and coefficients. Most useful for evaluating functions interior to finite elements

INPUT:
  coef:    [Vector{Float64}]
              Coefficients for some quantity evaluated at some point
  basis:   [Vector{Float64}]
              Nodal Basis functions evaluated at some point

OUTPUT:
  interp:  [Float64]
              Interpolated value
"""
function shapeEval(coef::Vector{T},basis::Vector{T}) where T<:Real
  interp  = zero(T)

  for i=1:length(basis)
    interp += coef[i]*basis[i]
  end

  return interp
end

"""
  derivShape2D(s,t,xNodes,yNodes,order)

Evaluates derivatives of shape functions

INPUT:
  s:       [Float64]
              Horizontal component of point in computational domain [-1,1]^2
  t:       [Float64]
              Vertical component of point in computational domain [-1,1]^2
  xNodes:  [Vector{Float64}]
              Horizontal component of element nodes
  yNodes:  [Vector{Float64}]
              Vertical component of element nodes
  order:   [Int]
              Order of basis polynomials on element

OUTPUT:
  phi:     [Vector{Float64}]
              Nodal basis function evaluated at node (s,t)
  dphidx:  [Vector{Float64}]
              Derivative in x-direction of nodal basis function
                evaluated at (s,t)
  dphidy:  [Vector{Float64}]
              Derivative in y-direction of nodal basis function
                evaluated at (s,t)
  jac:     [Float64]
              Jacobian of transformation from computation to physical
                domain (or is it the other way around...?)
"""
function derivShape2D(s::Float64,t::Float64,xNodes,yNodes,order::Int)
  if order==1
    numNodes=4
  elseif order==2
    numNodes=9
  else
    error("invalid order. Must be 1 or 2")
  end

  # evaluate shape functions
  phi,dphids,dphidt = shape2D(s,t,order)

  dxds = shapeEval(xNodes,dphids)
  dxdt = shapeEval(xNodes,dphidt)
  dyds = shapeEval(yNodes,dphids)
  dydt = shapeEval(yNodes,dphidt)

  # evaluates jacobian
  jac = dxds*dydt - dxdt*dyds

  # check element jacobian
  if jac < 1e-9
    println("Bad element warning ...\n")

    if jac <= 0.0
      error("singular Jacobian ... Aborted ...")
    end
  end

  # evaluate derivatives
  dphidx = @.  dphids*dydt - dphidt*dyds
  dphidy = @. -dphids*dxdt + dphidt*dxds

  return phi,dphidx,dphidy,jac
end

"""
  Evaluates dphidx and dphidy by reference but ALSO return jacobian float
"""
function derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s,t,xNodes,yNodes,order::Int)
(order==1 ? numNodes=4 : numNodes=9)

  dxds = shapeEval(xNodes,dphids)
  dxdt = shapeEval(xNodes,dphidt)
  dyds = shapeEval(yNodes,dphids)
  dydt = shapeEval(yNodes,dphidt)

  # evaluate derivatives
  for i=1:numNodes
  dphidx[i] =  dphids[i]*dydt - dphidt[i]*dyds
  dphidy[i] = -dphids[i]*dxdt + dphidt[i]*dxds
  end

  nothing
end

function jacCalc(xNodes,yNodes,dphids,dphidt)
  dxds = shapeEval(xNodes,dphids)
  dxdt = shapeEval(xNodes,dphidt)
  dyds = shapeEval(yNodes,dphids)
  dydt = shapeEval(yNodes,dphidt)

  # evaluates jacobian
  jac = dxds*dydt - dxdt*dyds

  # check element jacobian
  if jac < 1e-9
    println("Bad element warning ...\n")

    if jac <= 0.0
      error("singular Jacobian ... Aborted ...")
    end
  end

  return jac
end

"""
  ~~~  EFFICIENT ~~~

  GaussQuadPoints2D(order::Int)

Produces Gauss Quadrature weights and points, normalized to (-1,1), for a given order

INPUT:
  order:  [Int]
             Order of gauss quadrature

OUTPUT:
  w:      [Vector{Float64}]
             Weights of Guass Quadrature values
  s:      [Vector{Float64}]
             Horizontal components of Gauss Quadrature nodes
  t:      [Vector{Float64}]
             Horizontal components of Gauss Quadrature nodes
"""
function GaussQuadPoints2D(order::Int)
  if order==2        # 2x2 Gauss points
    w = Array{Float64,1}(4)
    s = Array{Float64,1}(4)
    t = Array{Float64,1}(4)

    gpt=1.0e0/sqrt(3.0e0);
    s[1] = -gpt; t[1] = -gpt; w[1]=1.0;
    s[2] =  gpt; t[2] = -gpt; w[2]=1.0;
    s[3] =  gpt; t[3] =  gpt; w[3]=1.0;
    s[4] = -gpt; t[4] =  gpt; w[4]=1.0;
  elseif order==3    # 3x3 Gauss points
    w = zeros(Float64,9)
    s = zeros(Float64,9)
    t = zeros(Float64,9)

    gpt=sqrt(0.6);
    s[1] = -gpt; t[1] = -gpt; w[1]=25/81
    s[2] =  gpt; t[2] = -gpt; w[2]=25/81
    s[3] =  gpt; t[3] =  gpt; w[3]=25/81
    s[4] = -gpt; t[4] =  gpt; w[4]=25/81
    s[5] =  0.0; t[5] = -gpt; w[5]=40/81
    s[6] =  gpt; t[6] =  0.0; w[6]=40/81
    s[7] =  0.0; t[7] =  gpt; w[7]=40/81
    s[8] = -gpt; t[8] =  0.0; w[8]=40/81
    s[9] =  0.0; t[9] =  0.0; w[9]=64/81
  elseif order==4      # 4x4 Gauss points
    w = zeros(Float64,16)
    s = zeros(Float64,16)
    t = zeros(Float64,16)

    tmp  = 2/7*sqrt(6/5)
    gpt1 = sqrt(3/7-tmp)
    gpt2 = sqrt(3/7+tmp)
    wp   = (18+sqrt(30))/36
    wm   = (18-sqrt(30))/36
    s[1]  =  gpt1; t[1]  =  gpt1; w[1]  = wp*wp
    s[2]  =  gpt1; t[2]  = -gpt1; w[2]  = wp*wp
    s[3]  =  gpt1; t[3]  =  gpt2; w[3]  = wp*wm
    s[4]  =  gpt1; t[4]  = -gpt2; w[4]  = wp*wm
    s[5]  = -gpt1; t[5]  =  gpt1; w[5]  = wp*wp
    s[6]  = -gpt1; t[6]  = -gpt1; w[6]  = wp*wp
    s[7]  = -gpt1; t[7]  =  gpt2; w[7]  = wp*wm
    s[8]  = -gpt1; t[8]  = -gpt2; w[8]  = wp*wm

    s[9]  =  gpt2; t[9]  =  gpt1; w[9]  = wm*wp
    s[10] =  gpt2; t[10] = -gpt1; w[10] = wm*wp
    s[11] =  gpt2; t[11] =  gpt2; w[11] = wm*wm
    s[12] =  gpt2; t[12] = -gpt2; w[12] = wm*wm
    s[13] = -gpt2; t[13] =  gpt1; w[13] = wm*wp
    s[14] = -gpt2; t[14] = -gpt1; w[14] = wm*wp
    s[15] = -gpt2; t[15] =  gpt2; w[15] = wm*wm
    s[16] = -gpt2; t[16] = -gpt2; w[16] = wm*wm
  else
    error("Check Gauss point integration specification")
  end

  return w,s,t
end
function GaussQuadPoints2DBigFloat(order::Int)
  if order==2        # 2x2 Gauss points
    w = Array{BigFloat,1}(4)
    s = Array{BigFloat,1}(4)
    t = Array{BigFloat,1}(4)

    gpt=1.0e0/sqrt(3.0e0);
    s[1] = -gpt; t[1] = -gpt; w[1]=1.0;
    s[2] =  gpt; t[2] = -gpt; w[2]=1.0;
    s[3] =  gpt; t[3] =  gpt; w[3]=1.0;
    s[4] = -gpt; t[4] =  gpt; w[4]=1.0;
  elseif order==3    # 3x3 Gauss points
    w = zeros(BigFloat,9)
    s = zeros(BigFloat,9)
    t = zeros(BigFloat,9)

    gpt=BigFloat(sqrt(0.6))
    s[1] = -gpt;       t[1] = -gpt;        w[1]=BigFloat(25/81)
    s[2] =  gpt;       t[2] = -gpt;        w[2]=BigFloat(25/81)
    s[3] =  gpt;       t[3] =  gpt;        w[3]=BigFloat(25/81)
    s[4] = -gpt;       t[4] =  gpt;        w[4]=BigFloat(25/81)
    s[5] =  zero(gpt); t[5] = -gpt;        w[5]=BigFloat(40/81)
    s[6] =  gpt;       t[6] =  zero(gpt);  w[6]=BigFloat(40/81)
    s[7] =  zero(gpt); t[7] =  gpt;        w[7]=BigFloat(40/81)
    s[8] = -gpt;       t[8] =  zero(gpt);  w[8]=BigFloat(40/81)
    s[9] =  zero(gpt); t[9] =  zero(gpt);  w[9]=BigFloat(64/81)
  elseif order==4      # 4x4 Gauss points
    w = zeros(BigFloat,16)
    s = zeros(BigFloat,16)
    t = zeros(BigFloat,16)

    tmp  = 2/7*sqrt(6/5)
    gpt1 = sqrt(3/7-tmp)
    gpt2 = sqrt(3/7+tmp)
    wp   = (18+sqrt(30))/36
    wm   = (18-sqrt(30))/36
    s[1]  =  gpt1; t[1]  =  gpt1; w[1]  = wp*wp
    s[2]  =  gpt1; t[2]  = -gpt1; w[2]  = wp*wp
    s[3]  =  gpt1; t[3]  =  gpt2; w[3]  = wp*wm
    s[4]  =  gpt1; t[4]  = -gpt2; w[4]  = wp*wm
    s[5]  = -gpt1; t[5]  =  gpt1; w[5]  = wp*wp
    s[6]  = -gpt1; t[6]  = -gpt1; w[6]  = wp*wp
    s[7]  = -gpt1; t[7]  =  gpt2; w[7]  = wp*wm
    s[8]  = -gpt1; t[8]  = -gpt2; w[8]  = wp*wm

    s[9]  =  gpt2; t[9]  =  gpt1; w[9]  = wm*wp
    s[10] =  gpt2; t[10] = -gpt1; w[10] = wm*wp
    s[11] =  gpt2; t[11] =  gpt2; w[11] = wm*wm
    s[12] =  gpt2; t[12] = -gpt2; w[12] = wm*wm
    s[13] = -gpt2; t[13] =  gpt1; w[13] = wm*wp
    s[14] = -gpt2; t[14] = -gpt1; w[14] = wm*wp
    s[15] = -gpt2; t[15] =  gpt2; w[15] = wm*wm
    s[16] = -gpt2; t[16] = -gpt2; w[16] = wm*wm
  else
    error("Check Gauss point integration specification")
  end

  return w,s,t

end

"""
  GaussEdge(edgeNumber,order)

Returns weights and point values for edges of gauss points
edgeNumber = 1 -> between points 1 and 2
edgeNumber = 2 -> between points 2 and 3
edgeNumber = 3 -> between points 3 and 4
edgeNumber = 4 -> between points 4 and 1

INPUT:
  edgeNumber:  [Int]
                  Index to determine where edge of interest lies
  order:       [Int]
                  Order of basis interpolating function

OUTPUT:
  w:           [Vector{Float64}]
                  Weights of Guass Quadrature values
  s:           [Vector{Float64}]
                  Horizontal components of Gauss Quadrature nodes
  t:           [Vector{Float64}]
                  Horizontal components of Gauss Quadrature nodes
"""
function GaussEdge(edgeNumber::Int,order::Int)
  w,_ = GaussQuadPoints1D(order)

  if edgeNumber == 1
    _,s = GaussQuadPoints1D(order)
    t   = -ones(order+1)
  elseif edgeNumber==2
    s   = ones(order+1)
    _,t = GaussQuadPoints1D(order)
  elseif edgeNumber==3
    _,s = GaussQuadPoints1D(order)
    t   = ones(order+1)
  elseif edgeNumber==4
    s   = -ones(order+1)
    _,t = GaussQuadPoints1D(order)
  else
    error("Invalid Edge Number entered. Must be 1-4, but $(edgeNumber) was entered")
  end

  return w::Vector{Float64},s::Vector{Float64},t::Vector{Float64}
end


"""
  FEMstats(xy,cm,dNodes,nNodes,F)

Prints to the terminal useful information about the system being solved

INPUT:
  mesh:   [AbstractMesh]
             Mesh object containing nodes and connectivity
  prob:   [AbstractProblem]
             Problem object containing boundary information
  LinOp:  [AbstractLinearOperator]
             Linear Operator object containing discretized operator matrix and
               forcing vector

OUTPUT:
  none
"""
function FEMstats(mesh,prob,LinOp)
  numTnodes = length(mesh.xy)
  Nel       = length(mesh.cm)
  numDnodes = length(prob.bcid[:dNodes])
  numNnodes = length(prob.bcid[:nNodes])
  sysSize   = length(LinOp.F)

  println("# of total nodes:        ",numTnodes)
  println("# of elements:           ",Nel)
  println("# of Dirichlet BC nodes: ",numDnodes)
  println("# of Neumann BC nodes:   ",numNnodes)
  println("% boundary nodes:        ",
          round(numTnodes/(numDnodes+numNnodes),2),"%")
  println("Size of linear system:   ",sysSize," x ",sysSize)
end

"""
	GaussQuadPoints1D(order)

Produces Gauss Quadrature weights and points, normalized to (-1,1), for a given order
"""
function GaussQuadPoints1D(order::Int)
	w = zeros(Float64,order)
	s = zeros(Float64,order)

	if order == 2
		w[1] = 1
		w[2] = 1
		s[1] = sqrt(1/3)
		s[2] = -s[1]
  elseif order == 3
    w[1] = 8/9
    w[2] = 5/9
    w[3] = w[2]
    s[1] = 0.0
    s[2] =  sqrt(3/5)
    s[3] = -s[2]
  end

  return w, s
end
