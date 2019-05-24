# File loaded into eFEMpart

"""
  shape2D(s::Real,t::Real,order::Int)

Computes shape functions of linear or quadratic order on a quadrilateral element

ASSUMPTIONS:
- Only available for quadratic elements
- Assumes computing with double-precision arithmetic

INPUT:
  s:      [Real]
              horizontal location in computation domain [-1,1]^2
  t:      [Real]
              vertical location in computation domain [-1,1]^2
  order:  [Int]
              Order of basis function (1 -> linear, 2-> quadratic)

OUTPUT:
  phi:    [Vector{Real}]
              Nodal basis functions evaluated at (s,t)
  dphids  [Vector{Real}]
              Derivative in s-direction of nodal basis function
                evaluated at (s,t)
  dphidt:  [Vector{Real}]
              Derivative in t-direction of nodal basis function
                evaluated at (s,t)
"""
function shape2D(s::T,t::T,order::Int) where T<:Real
  ein = one(s)
  tt = typeof(s)

  # linear element
  if order == 1
    phi    = Array{tt}(undef,4)
    dphids = Array{tt}(undef,4)
    dphidt = Array{tt}(undef,4)

    phi[1]    = 0.25*(ein-s)*(ein-t)
    phi[2]    = 0.25*(ein+s)*(ein-t)
    phi[3]    = 0.25*(ein+s)*(ein+t)
    phi[4]    = 0.25*(ein-s)*(ein+t)
    dphids[1] = -0.25*(ein-t)
    dphids[2] =  0.25*(ein-t)
    dphids[3] =  0.25*(ein+t)
    dphids[4] = -0.25*(ein+t)
    dphidt[1] = -0.25*(ein-s)
    dphidt[2] = -0.25*(ein+s)
    dphidt[3] =  0.25*(ein+s)
    dphidt[4] =  0.25*(ein-s)
  
  # quadratic element
  elseif order == 2
    phi    = Array{tt}(undef,9)
    dphids = Array{tt}(undef,9)
    dphidt = Array{tt}(undef,9)

    ells1  = 0.5*s*(s-ein); ellt1 = 0.5*t*(t-ein);
    ells2  = ein-(s*s);     ellt2 = ein-(t*t);
    ells3  = 0.5*s*(s+ein); ellt3 = 0.5*t*(t+ein);
    dells1 = s-0.5;         dellt1 = t-0.5;
    dells2 = -2*s;          dellt2 = -2*t;
    dells3 = s+0.5;         dellt3 = t+0.5;

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
  shape2D!(phi::Vector{T},dphids::Vector{T},dphidt::Vector{T},
           s::T,t::T,order::Int) where T<:Real

in-place version of shape2D; see shape2D function info.
"""
function shape2D!(phi,dphids,dphidt,s::T,t::T,order::Int) where T<:Real
  ein = one(s)

  # linear element
  if order == 1
    phi[1]    = 0.25*(ein-s)*(ein-t)
    phi[2]    = 0.25*(ein+s)*(ein-t)
    phi[3]    = 0.25*(ein+s)*(ein+t)
    phi[4]    = 0.25*(ein-s)*(ein+t)
    dphids[1] = -0.25*(ein-t)
    dphids[2] =  0.25*(ein-t)
    dphids[3] =  0.25*(ein+t)
    dphids[4] = -0.25*(ein+t)
    dphidt[1] = -0.25*(ein-s)
    dphidt[2] = -0.25*(ein+s)
    dphidt[3] =  0.25*(ein+s)
    dphidt[4] =  0.25*(ein-s)

  # quadratic element
  elseif order == 2
    ells1  = 0.5*s*(s-ein);  ellt1 = 0.5*t*(t-ein);
    ells2  = ein-(s*s);      ellt2 = ein-(t*t);
    ells3  = 0.5*s*(s+ein);  ellt3 = 0.5*t*(t+ein);
    dells1 = s-0.5;          dellt1 = t-0.5;
    dells2 = -2*s;           dellt2 = -2*t;
    dells3 = s+0.5;          dellt3 = t+0.5;

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
  shapeEval(coef::Vector{Real},basis::Vector{Real})

Evaluates the sum given by the basis values and coefficients. Most useful for 
evaluating functions interior to finite elements when (s,t) are known

INPUT:
  coef:    [Vector{Real}]
              Coefficients for some quantity evaluated at some point
  basis:   [Vector{Real}]
              Nodal Basis functions evaluated at some point

OUTPUT:
  -        [Real]
              Interpolated value
"""
function shapeEval(coef::Vector{T},basis::Vector{T}) where T<:Real
  return LinearAlgebra.dot(coef,basis)
end

"""
  derivShape2D(s,t,xNodes,yNodes,order)

Evaluates derivatives of shape functions

INPUT:
  s:       [Real]
              Horizontal component of point in computational domain [-1,1]^2
  t:       [Real]
              Vertical component of point in computational domain [-1,1]^2
  xNodes:  [Vector{Real}]
              Horizontal component of element nodes
  yNodes:  [Vector{Real}]
              Vertical component of element nodes
  order:   [Int]
              Order of basis polynomials on element

OUTPUT:
  phi:     [Vector{Real}]
              Nodal basis function evaluated at node (s,t)
  dphidx:  [Vector{Real}]
              Derivative in x-direction of nodal basis function
                evaluated at (s,t)
  dphidy:  [Vector{Real}]
              Derivative in y-direction of nodal basis function
                evaluated at (s,t)
  jac:     [Real]
              Jacobian of transformation from computation to physical
                domain (or is it the other way around...?)
"""
function derivShape2D(s::T,t::T,xNodes::Vector{T},yNodes::Vector{T},order::Int) where T<:Real
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
  derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s,t,xNodes,yNodes,order)

Evaluates derivatives of shape functions in-place
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

"""
  jacCalc(xNodes,yNodes,dphids,dphidt)

Computes jacobian of element
"""
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
  GaussQuadPoints2D(order::Int)

Produces Gauss Quadrature weights and points, normalized to (-1,1), for a given order

INPUT:
  order:  [Int]
             Order of gauss quadrature

OUTPUT:
  w:      [Vector{Real}]
             Weights of Guass Quadrature values
  s:      [Vector{Real}]
             Horizontal components of Gauss Quadrature nodes
  t:      [Vector{Real}]
             Horizontal components of Gauss Quadrature nodes
"""
function GaussQuadPoints2D(order::Int;tt=Float64::DataType)
  if order==2        # 2x2 Gauss points
    w = Array{tt,1}(undef,4)
    s = Array{tt,1}(undef,4)
    t = Array{tt,1}(undef,4)

    gpt=1.0e0/sqrt(3.0e0);
    s[1] = -gpt; t[1] = -gpt; w[1]=1.0;
    s[2] =  gpt; t[2] = -gpt; w[2]=1.0;
    s[3] =  gpt; t[3] =  gpt; w[3]=1.0;
    s[4] = -gpt; t[4] =  gpt; w[4]=1.0;
  elseif order==3    # 3x3 Gauss points
    w = zeros(tt,9)
    s = zeros(tt,9)
    t = zeros(tt,9)

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
    w = zeros(tt,16)
    s = zeros(tt,16)
    t = zeros(tt,16)

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
  w:           [Vector{Real}]
                  Weights of Guass Quadrature values
  s:           [Vector{Real}]
                  Horizontal components of Gauss Quadrature nodes
  t:           [Vector{Real}]
                  Horizontal components of Gauss Quadrature nodes
"""
function GaussEdge(edgeNumber::Int,order::Int;tt=Float64::DataType)
  w,_ = GaussQuadPoints1D(order;tt=tt)

  if edgeNumber == 1
    _,s = GaussQuadPoints1D(order;tt=tt)
    t   = -ones(order+1)
  elseif edgeNumber==2
    s   = ones(order+1)
    _,t = GaussQuadPoints1D(order;tt=tt)
  elseif edgeNumber==3
    _,s = GaussQuadPoints1D(order;tt=tt)
    t   = ones(order+1)
  elseif edgeNumber==4
    s   = -ones(order+1)
    _,t = GaussQuadPoints1D(order;tt=tt)
  else
    error("Invalid Edge Number entered. Must be 1-4, but $(edgeNumber) was entered")
  end

  return w::Vector{tt},s::Vector{tt},t::Vector{tt}
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

function FEMstats(mesh)
  numTnodes = length(mesh.xy)
  Nel       = length(mesh.cm)

  println("# of total nodes:        ",numTnodes)
  println("# of elements:           ",Nel)
end

"""
	GaussQuadPoints1D(order)

Produces Gauss Quadrature weights and points, normalized to (-1,1), for a given order

Not really sure where this function is used...
"""
function GaussQuadPoints1D(order::Int;tt=Float64::DataType)
	w = zeros(tt,order)
	s = zeros(tt,order)

	if order == 2
		w[1] = one(tt)
		w[2] = one(tt)
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
