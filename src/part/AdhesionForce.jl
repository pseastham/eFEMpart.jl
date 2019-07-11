# Functions for computing particle-wall contact forces of 
# via quadrature of lennard-jones potential forces

include("CellLists.jl")      # loads in Point2D type
include("ParticleTypes.jl")  # load Wall types
include("CohesionForce.jl")  # loads in LJ force function
include("isInside.jl")       # loads onSegment()

### Compute nearest point to particle on LineWall, CircleWall, and ArcWall
function NearestPoint!(point::Point2D,node::Point2D,wall::LineWall)
    px=node.x; py=node.y

    Ax=wall.nodes[1].x; Ay=wall.nodes[1].y
    Bx=wall.nodes[2].x; By=wall.nodes[2].y

    bx=px-Ax; by=py-Ay
    ax=Bx-Ax; ay=By-Ay

    ℓ2 = ax^2+ay^2

    dotprod = ax*bx + ay*by
    
    point.x = dotprod*ax/ℓ2 + Ax
    point.y = dotprod*ay/ℓ2 + Ay

    nothing
end
function NearestPoint!(point::Point2D,node::Point2D,wall::CircleWall)
    px=node.x; py=node.y
    cx=wall.center.x; cy=wall.center.y
    r = wall.radius
    θ = atan(py-cy,px-cx)

    point.x = cx + r*cos(θ)
    point.y = cy + r*sin(θ)

    nothing
end
function NearestPoint!(point::Point2D,node::Point2D,wall::ArcWall)
    px= node.x; py=node.y
    cx= wall.nodes[2].x; cy=wall.nodes[2].y
    r = sqrt((cx-wall.nodes[1].x)^2 + (cy-wall.nodes[1].y)^2)
    θ = atan(py-cy,px-cx)

    point.x = cx + r*cos(θ)
    point.y = cy + r*sin(θ)

    nothing
end
function NearestPoint(node::Point2D,wall::W) where W<:AbstractWall
    point = Point2D(0.0,0.0)
    NearestPoint!(point,node,wall)
    return point
end

# checks whether node is close-enough to point (used to see if calculations necessary)
# radius is radius of particle 'node'
# k: number of times away from which quadrature should be considered
function isCloseEnough(point::Point2D,node::Point2D,k::T,rc::T,radius::T) where T<:Real
    ℓ = sqrt((point.x-node.x)^2 + (point.y-node.y)^2)

    return  (ℓ < k*rc*2*radius ? true : false)
end

# generates quadrature nodes for linewall, circlewall, and arcwall
# point: center point about which quadrature nodes are taken
# k: number of times away from which quadrature should be considered
# N: number of quadrature nodes
# xquad,yquad: quadrature nodes to be passed in to be computed
function generateQuadNodes!(xquad::Vector{T},yquad::Vector{T},point::Point2D,wall::LineWall,k::T,particleRadius::T) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end

    N = length(xquad)

    tx=wall.t[1]; ty=wall.t[2]
    L = 2*particleRadius
    Δs = 2*L/(N-1)

    temp1X = L*tx; temp2X = Δs*tx
    temp1Y = L*ty; temp2Y = Δs*ty

    for ti=1:N
        xquad[ti] = point.x - temp1X + (ti-1)*temp2X
        yquad[ti] = point.y - temp1Y + (ti-1)*temp2Y
    end

    nothing
end
function generateQuadNodes!(xquad::Vector{T},yquad::Vector{T},point::Point2D,wall::CircleWall,k::T,particleRadius::T) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end

    N = length(xquad)

    px=point.x; py=point.y
    cx=wall.center.x; cy=wall.center.y
    r = wall.radius
    θcenter = atan(py-cy,px-cx)
    θmin = θcenter - pi/2
    Δθ = pi/(N-1)

    for ti=1:N
        xquad[ti] = cx + r*cos(θmin+(ti-1)*Δθ)
        yquad[ti] = cy + r*sin(θmin+(ti-1)*Δθ)
    end

    nothing
end
function generateQuadNodes!(xquad::Vector{T},yquad::Vector{T},point::Point2D,wall::ArcWall,k::T,particleRadius::T) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end

    N = length(xquad)

    px=point.x; py=point.y
    cx= wall.nodes[2].x; cy=wall.nodes[2].y
    r = sqrt((cx-wall.nodes[1].x)^2 + (cy-wall.nodes[1].y)^2)
    θcenter = atan(py-cy,px-cx)
    θmin = θcenter - pi/2
    Δθ = 2*pi/(N-1)

    for ti=1:N
        xquad[ti] = cx + r*cos(θmin+(ti-1)*Δθ)
        yquad[ti] = cy + r*sin(θmin+(ti-1)*Δθ)
    end

    nothing
end

# function to determine whether quadrature node (sx,sy) is within line
function isInLine(wall::LineWall,sx::T,sy::T,s::Point2D) where T<:Real 
    s.x=sx; s.y=sy
    return onSegment(wall.nodes[1],s,wall.nodes[2])                         # onSegment is located in isInside.jl
end
# note: s is input only to make all arguments for isInLine the same
function isInLine(wall::CircleWall,sx::T,sy::T,s::Point2D) where T<:Real
    cx = wall.center.x; cy=wall.center.y
    val = abs(wall.radius - sqrt((cx-sx)^2 + (cy-sy)^2))
    TOL = 1e-12
    return (val < TOL ? true : false)
end
# note: s is input only to make all arguments for isInLine the same
function isInLine(wall::ArcWall,sx::T,sy::T,s::Point2D) where T<:Real
    cx = wall.nodes[2].x; cy=wall.nodes[2].y
    radius = sqrt((cx-wall.nodes[1].x)^2 + (cy-wall.nodes[1].y)^2)
    val = abs(radius - sqrt((cx-sx)^2 + (cy-sy)^2))

    θ1 = atan(wall.nodes[1].y-cy,wall.nodes[1].x-cx)
    θ2 = atan(wall.nodes[3].y-cy,wall.nodes[3].x-cx)
    θs = atan(sy-cy,sx-cx)

    TOL = 1e-12
    isOn = false
    if val < TOL
        if θ1 > θ2
            if θs > θ1 || θs < θ2
                isOn = true
            end
        else
            if θs > θ1 && θs < θ2
                return true
            end
        end
    end

    return isOn
end

# computes quadrature on wall
# node = position of particle that you are comparing
# ϵ = strength of lennard jones potential
# radius = particle radius
# xquad,yquad = quadrature points
# rc = cutoff percentage 
function wallTrapQuad(node::Point2D,wall::LineWall,
                      xquad::Vector{T},yquad::Vector{T},
                      ϵ::T,radius::T,rc::T,s::Point2D) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end    

    N = length(xquad)
    Δs = sqrt((wall.nodes[1].x - wall.nodes[2].x)^2 + (wall.nodes[1].y - wall.nodes[2].y)^2)/(N-1)
    
    # this should change to radius + 0.5*wall.width once this feature
    # is added to walls datatype
    d = 2*radius

    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=2:N
        if isInLine(wall,xquad[ti],yquad[ti],s)
            # compute vector from wall to node 
            qxKm1 = node.x-xquad[ti-1]; qyKm1 = node.y-yquad[ti-1]
            qxK   = node.x-xquad[ti];   qyK   = node.y-yquad[ti]

            # compute strength of LJ force
            FxKm1,FyKm1 = ForceCalculation(ϵ,d,rc,qxKm1,qyKm1)
            FxK,FyK     = ForceCalculation(ϵ,d,rc,qxK,qyK)

            quadSumX += 0.5*(FxKm1 + FxK)*Δs
            quadSumY += 0.5*(FyKm1 + FyK)*Δs
        end
    end

    return quadSumX,quadSumY
end
function wallTrapQuad(node::Point2D,wall::CircleWall,
                        xquad::Vector{T},yquad::Vector{T},
                        ϵ::T,particleradius::T,rc::T,s::Point2D) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end    

    N = length(xquad)
    cx = wall.center.x; cy = wall.center.y
    θ1 = atan(yquad[1]-cy,xquad[1]-cx)
    θ2 = atan(yquad[2]-cy,xquad[2]-cx)
    Δs = abs(wall.radius*(θ2-θ1))
    # this should change to radius + 0.5*wall.width once this feature
    # is added to walls datatype
    d = 2*particleradius

    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=1:N
        if isInLine(wall,xquad[ti],yquad[ti],s)
            # compute vector from x,y to p 
            qx=node.x - xquad[ti]; qy=node.y - yquad[ti]

            # compute strength of LJ force
            Fx,Fy = ForceCalculation(ϵ,d,rc,qx,qy)

            # add to running quadrature sum
            if (ti==1 || ti==N)
                quadSumX += 0.5*Fx
                quadSumY += 0.5*Fy
            else
                quadSumX += Fx
                quadSumY += Fy
            end
        end
    end

    quadSumX *= Δs
    quadSumY *= Δs

    return quadSumX,quadSumY
end
function wallTrapQuad(node::Point2D,wall::ArcWall,
                        xquad::Vector{T},yquad::Vector{T},
                        ϵ::T,particleradius::T,rc::T,s::Point2D) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end    

    N = length(xquad)
    cx = wall.nodes[2].x; cy = wall.nodes[2].y
    θ1 = atan(yquad[1]-cy,xquad[1]-cx)
    θ2 = atan(yquad[2]-cy,xquad[2]-cx)
    radius = sqrt((cx - wall.nodes[1].x)^2 + (cy - wall.nodes[1].y)^2)
    # this should change to radius + 0.5*wall.width once this feature
    # is added to walls datatype
    d = 2*particleradius

    if θ1 > θ2
        θ2 += 2pi
    end

    Δs = abs(radius*(θ2-θ1))
    N = length(xquad)
    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=1:N
        if isInLine(wall,xquad[ti],yquad[ti],s)
            # compute vector from x,y to p 
            qx=node.x - xquad[ti]; qy=node.y - yquad[ti]

            # compute strength of LJ force
            Fx,Fy = ForceCalculation(ϵ,d,rc,qx,qy)

            # add to running quadrature sum
            if (ti==1 || ti==N)
                quadSumX += 0.5*Fx
                quadSumY += 0.5*Fy
            else
                quadSumX += Fx
                quadSumY += Fy
            end
        end
    end

    quadSumX *= Δs
    quadSumY *= Δs

    return quadSumX,quadSumY
end

# Computes adhesion force for all particles against all possible walls
function AdhesionForce!(afX::Vector{T},afY::Vector{T},pList::Vector{Point2D},rList::Vector{T},wList::Vector{W},
                        k::T,rc::T,ϵ::T,pointOnWall::Point2D,xquad::Vector{T},yquad::Vector{T}) where {W<:AbstractWall, T<:Real}
    if length(afX) != length(afY)
        throw(DimensionMismatch("length of afX and afY must match!"))
    end

    for tp=1:length(pList), tw=1:length(wList)
        NearestPoint!(pointOnWall,pList[tp],wList[tw])

        if isCloseEnough(pointOnWall,pList[tp],k,rc,rList[tp])
            generateQuadNodes!(xquad,yquad,pointOnWall,wList[tw],k,rList[tp])
            Fx,Fy = wallTrapQuad(pList[tp],wList[tw],xquad,yquad,ϵ,rList[tp],rc,pointOnWall)
            afX[tp] += Fx
            afY[tp] += Fy
        end 
    end

    nothing
end