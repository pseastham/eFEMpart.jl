# File to be loaded into eFEMpart

# -- Types for surface quadratures
# can probably be summarized somewhere else, esp Point2D
# used in "part" component of code
mutable struct Point 
    x::Float64
    y::Float64
end

mutable struct mySurface
    P::Point
    Q::Point
end

function mySurface(Parr::Vector{Float64},Qarr::Vector{Float64})
  P = Point(Parr[1],Parr[2])
  Q = Point(Qarr[1],Qarr[2])
  return mySurface(P,Q)
end

"""
  DomainQuad(mesh,u)

  Computes the 2D integral over the entire domain defined in mesh of the 
  solution u. Uses Gaussian Quadrature.
"""
function DomainQuad(mesh::M,u::Vector{Float64}) where {M<:AbstractMesh}
  Nel   = length(mesh.cm)
  ngpt  = length(mesh.cm[1].NodeList)
  order = (ngpt==4 ? 1 : 2)
  quad  = 0.0
  
  # generate basis values at gauss points
  w,s,t = GaussQuadPoints2D(order+1)

  xN     = zeros(Float64,ngpt)
  yN     = zeros(Float64,ngpt)
  ucoefs = zeros(Float64,ngpt)

  for el=1:Nel
    for ti=1:ngpt
      c          = mesh.cm[el].NodeList[ti]
      xN[ti]     = mesh.xy[c].x
      yN[ti]     = mesh.xy[c].y
      ucoefs[ti] = u[c]
    end
    
    for gpt=1:ngpt
      phi,_,_,jac = derivShape2D(s[gpt],t[gpt],xN,yN,order)
      ugpt = shapeEval(ucoefs,phi)
      
      quad += w[gpt]*abs(ugpt)*jac
    end
  end
  
  return quad
end
function DomainQuad(mesh,u::S) where S<:AbstractSolution
  uarr = u.u
  return DomainQuad(mesh,uarr)
end

"""
  SurfaceQuad(mesh,u,surf,N)

  Computes N-discretized line integral of u on 'surface' (aka line) surf. 
  Quadrature uses trapezoid rule 
"""
function SurfaceQuad(mesh::M,u::Vector{Float64},surf::mySurface,N::Int) where 
                    {M<:AbstractMesh}
  (mesh.order==:Linear ? order=1 : order=2)
  (order==1 ? nGpt=2 : nGpt=3)

  # Define 2 points
  P = surf.P; Q = surf.Q
  rlen = sqrt((P.x-Q.x)^2 + (P.y-Q.y)^2)

  # Interpolate N segments (N+1 points)
  Δt = rlen/N
  xN = zeros(Float64,N+1)
  yN = zeros(Float64,N+1)
  uN = zeros(Float64,N+1)
  tx,ty = ComputeTangent(P,Q)

  xN[1] = P.x; yN[1] = P.y
  for i=2:(N+1)
    xN[i] = xN[1] + (i-1)*Δt*tx
    yN[i] = yN[1] + (i-1)*Δt*ty
  end

  # Interpolate solution at N+1 points
  for i=1:(N+1)
    point = [xN[i];yN[i]]
    uN[i] = pointTransform(mesh,u,point)
  end

  # For each segment, do trapezoid rule quadrature
  quad = 0.0
  for ti=2:N
    quad += uN[ti]*Δt
  end
  quad += 0.5*Δt*(uN[1] + uN[N+1])

  # return quad
  return quad
end
function SurfaceQuad(mesh::M,u::S,surf::mySurface,N::Int) where 
  {M<:AbstractMesh,S<:AbstractSolution}
  uarr = u.u
  return SurfaceQuad(mesh,uarr,surf,N) 
end

"""
  for some scalar field c, computes grad(c)*n across a surface
"""
function SurfaceFlux(mesh::M,u::Vector{Float64},surf::mySurface,Nquad::Int) where 
                    {M<:AbstractMesh}
  Npts = length(mesh.xy)
  Nel  = length(mesh.cm)
  (mesh.order==:Linear ? order=1 : order=2)
  (order==1 ? nGpt=2 : nGpt=3)
  NnodesPerElm = (order==1 ? 4 : 9)

  # Define 2 points
  P = surf.P; Q = surf.Q
  rlen = sqrt((P.x-Q.x)^2 + (P.y-Q.y)^2)

  # Initialize arrays
  Δt = rlen/Nquad
  xN = zeros(Float64,Nquad+1)
  yN = zeros(Float64,Nquad+1)
  uN = zeros(Float64,Nquad+1)

  dudx  = zeros(Float64,Npts)
  dudy  = zeros(Float64,Npts)
  dudxN = zeros(Float64,Nquad+1)
  dudyN = zeros(Float64,Nquad+1)

  tx,ty = ComputeTangent(P,Q)
  nx,ny = ComputeNormal(P,Q)

  # Interpolate N segments (N+1 points)
  xN[1] = P.x; yN[1] = P.y
  for i=2:(Nquad+1)
    xN[i] = xN[1] + (i-1)*Δt*tx
    yN[i] = yN[1] + (i-1)*Δt*ty
  end

  # compute derivatives at mesh points
  s = [-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,-1.0,0.0]
  t = [-1.0,-1.0,1.0,1.0,-1.0,0.0,1.0,0.0,0.0]
  for el=1:Nel
    cv = mesh.cm[el].NodeList
    xNMesh = [mesh.xy[i].x for i in cv]
    yNMesh = [mesh.xy[i].y for i in cv]
    utemp  = [u[i] for i in cv]

    for i=1:NnodesPerElm
      # compute basis at node points
      _,dphidx,dphidy,jac = derivShape2D(s[i],t[i],xNMesh,yNMesh,order)
      dudx[mesh.cm[el].NodeList[i]] = shapeEval(utemp,dphidx)/jac
      dudy[mesh.cm[el].NodeList[i]] = shapeEval(utemp,dphidy)/jac
    end
  end

  # Interpolate derivative solution at N+1 points
  for i=1:(Nquad+1)
    point  = [xN[i];yN[i]]
    dudxN[i]  = pointTransform(mesh,dudx,point)
    dudyN[i]  = pointTransform(mesh,dudy,point)
  end

  # For each segment, do trapezoid rule quadrature
  quad = 0.0
  for ti=2:Nquad
    quad += (dudxN[ti]*nx + dudyN[ti]*ny)*Δt
  end
  quad += 0.5*Δt*(dudxN[1]*nx + dudyN[1]*ny + dudxN[Nquad+1]*nx + dudyN[Nquad+1]*ny)

  # return quad
  return quad
end
function SurfaceFlux(mesh::M,u::S,surf::mySurface,N::Int) where 
  {M<:AbstractMesh,S<:AbstractSolution}
  uarr = u.u
  return SurfaceFlux(mesh,uarr,surf,N) 
end

function VelocityFlux(mesh::M,u::Vector{Float64},v::Vector{Float64},
                surf::mySurface,N::Int) where {M<:AbstractMesh}
  # compute normal
  nx,ny = ComputeNormal(surf.P,surf.Q)

  # copmute u-flux
  uFlux = SurfaceQuad(mesh,u,surf,N)

  # scale u-flux by horizontal normal
  uval = nx*uFlux

  # compute v-flux
  vFlux = SurfaceQuad(mesh,v,surf,N)

  # scale v-flux by vertical normal
  vval = ny*vFlux

  return uval + vval
end

#############################################
# NECESSARY FUNCTIONS FOR SURFACE Quadrature
#############################################

"""
  Interp2Dline(surf,N)

  Interpolates the points defining a surface in 2D (aka a line) into N segments.
"""
function Interp2Dline(surf::mySurface,N::Int)
  vt = ComputeTangent(surf.P,surf.Q)
  vn = ComputeNormal(surf.P,surf.Q)
  
  return vt,vn
end

# computes the normal at 2 points p1=[x1,y1] and p2=[x2,y2]
function ComputeTangent(p1,p2)
  rx = p2.x - p1.x
  ry = p2.y - p1.y

  rmag = sqrt(rx^2 + ry^2)

  tx = rx/rmag
  ty = ry/rmag

  return tx,ty
end

# computes the normal at 2 points p1=[x1,y1] and p2=[x2,y2]
# Normal is defined if you are standing on p1 and looking at p2, the normal goes to 
# the right
function ComputeNormal(p1,p2)
  tx,ty = ComputeTangent(p1,p2)

  nx = ty
  ny = -tx

  return nx,ny
end
