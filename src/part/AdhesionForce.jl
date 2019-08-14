# Functions for computing particle-wall contact forces of 
# via quadrature of lennard-jones potential forces

include("CellLists.jl")      # loads in Point2D type
include("ParticleTypes.jl")  # load Wall types
include("CohesionForce.jl")  # loads in LJ force function
include("isInside.jl")       # loads onSegment()

### Compute nearest point to particle on LineWall, CircleWall, and ArcWall
function NearestPoint!(point::Point2D{T},node::Point2D{T},wall::LineWall) where T<:Real
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
function NearestPoint!(point::Point2D{T},node::Point2D{T},wall::CircleWall) where T<:Real
    px=node.x; py=node.y
    cx=wall.center.x; cy=wall.center.y
    r = wall.radius
    θ = atan(py-cy,px-cx)

    point.x = cx + r*cos(θ)
    point.y = cy + r*sin(θ)

    nothing
end
function NearestPoint!(point::Point2D{T},node::Point2D{T},wall::ArcWall) where T<:Real
    px= node.x; py=node.y
    cx= wall.nodes[2].x; cy=wall.nodes[2].y
    r = sqrt((cx-wall.nodes[1].x)^2 + (cy-wall.nodes[1].y)^2)
    θ = atan(py-cy,px-cx)

    point.x = cx + r*cos(θ)
    point.y = cy + r*sin(θ)

    nothing
end
function NearestPoint(node::Point2D{T},wall::W) where {T<:Real,W<:AbstractWall}
    point = Point2D(0.0,0.0)
    NearestPoint!(point,node,wall)
    return point
end

# checks whether node is close-enough to point (used to see if calculations necessary)
# radius is radius of particle 'node'
# k: number of times away from which quadrature should be considered
function isCloseEnough(point::Point2D{T},node::Point2D{T},k,radius::T) where T<:Real
    ℓ = sqrt((point.x-node.x)^2 + (point.y-node.y)^2)

    return  (ℓ < k*2*radius ? true : false)
end

# generates quadrature nodes for linewall, circlewall, and arcwall
# point: center point about which quadrature nodes are taken
# k: number of times away from which quadrature should be considered
# N: number of quadrature nodes
# xquad,yquad: quadrature nodes to be passed in to be computed
function generateQuadNodes!(xquad::Vector{T},yquad::Vector{T},point::Point2D{T},node::Point2D{T},wall::LineWall{T},k,particleRadius::T) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end

    N = length(xquad)

    tx=wall.t[1]; ty=wall.t[2]
    dist2wall = sqrt((node.x - point.x)^2+(node.y-point.y)^2)
    L = 2*particleRadius/dist2wall
    Δs = 2*L/(N-1)

    temp1X = L*tx; temp2X = Δs*tx
    temp1Y = L*ty; temp2Y = Δs*ty

    for ti=1:N
        xquad[ti] = point.x - temp1X + (ti-1)*temp2X
        yquad[ti] = point.y - temp1Y + (ti-1)*temp2Y
    end

    nothing
end
function generateQuadNodes!(xquad::Vector{T},yquad::Vector{T},point::Point2D{T},node::Point2D{T},wall::CircleWall{T},k,particleRadius::T) where T<:Real
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
function generateQuadNodes!(xquad::Vector{T},yquad::Vector{T},point::Point2D{T},node::Point2D{T},wall::ArcWall{T},k,particleRadius::T) where T<:Real
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
function isInLine(wall::LineWall{T},sx::T,sy::T,s::Point2D{T}) where T<:Real 
    s.x=sx; s.y=sy
    return onSegment(wall.nodes[1],s,wall.nodes[2])                         # onSegment is located in isInside.jl
end
# note: s is input only to make all arguments for isInLine the same
function isInLine(wall::CircleWall{T},sx::T,sy::T,s::Point2D{T}) where T<:Real
    cx = wall.center.x; cy=wall.center.y
    val = abs(wall.radius - sqrt((cx-sx)^2 + (cy-sy)^2))
    TOL = 1e-12
    return (val < TOL ? true : false)
end
# note: s is input only to make all arguments for isInLine the same
function isInLine(wall::ArcWall{T},sx::T,sy::T,s::Point2D{T}) where T<:Real
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
# k = cutoff percentage 
function wallTrapQuad(node::Point2D{T},wall::LineWall{T},
                      xquad::Vector{T},yquad::Vector{T},
                      ϵ::T,particleradius::T,k,s::Point2D{T}) where T<:Real
    if length(xquad) != length(yquad)
        throw(DimensionMismatch("xquad and yquad must be same length"))
    end    

    N = length(xquad)
    Δs = sqrt((wall.nodes[1].x - wall.nodes[2].x)^2 + (wall.nodes[1].y - wall.nodes[2].y)^2)/(N-1)
    
    # this should change to radius + 0.5*wall.width once this feature
    # is added to walls datatype
    d = 1.0536*2*particleradius

    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=2:N
        if isInLine(wall,xquad[ti],yquad[ti],s)
            # compute vector from wall to node 
            qxKm1 = xquad[ti-1] - node.x; qyKm1 = yquad[ti-1] - node.y
            qxK   = xquad[ti]   - node.x;   qyK = yquad[ti]   - node.y

            # compute strength of LJ force
            FxKm1,FyKm1 = ForceCalculation(ϵ,d,k,qxKm1,qyKm1)
            FxK,FyK     = ForceCalculation(ϵ,d,k,qxK,qyK)

            quadSumX += 0.5*(FxKm1 + FxK)*Δs
            quadSumY += 0.5*(FyKm1 + FyK)*Δs
        end
    end

    return quadSumX,quadSumY
end
function wallTrapQuad(node::Point2D{T},wall::CircleWall{T},
                        xquad::Vector{T},yquad::Vector{T},
                        ϵ::T,particleradius::T,k,s::Point2D{T}) where T<:Real
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
    d = 1.05*2*particleradius

    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=1:N
        if isInLine(wall,xquad[ti],yquad[ti],s)
            # compute vector from x,y to p 
            qx=xquad[ti]-node.x; qy=yquad[ti]-node.y 

            # compute strength of LJ force
            Fx,Fy = ForceCalculation(ϵ,d,k,qx,qy)

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
function wallTrapQuad(node::Point2D{T},wall::ArcWall{T},
                        xquad::Vector{T},yquad::Vector{T},
                        ϵ::T,particleradius::T,k,s::Point2D{T}) where T<:Real
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
    d = 1.05*2*particleradius

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
            qx=xquad[ti]-node.x
            qy=yquad[ti]-node.y

            # compute strength of LJ force
            Fx,Fy = ForceCalculation(ϵ,d,k,qx,qy)

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
function AdhesionForce!(afX::Vector{T},afY::Vector{T},pList::Vector{Point2D{T}},rList::Vector{T},wList::Vector{W},
                        k,ϵ::T,pointOnWall::Point2D{T},xquad::Vector{T},yquad::Vector{T}) where {T<:Real,W<:AbstractWall}
    if length(afX) != length(afY)
        throw(DimensionMismatch("length of afX and afY must match!"))
    end

    # zero-out cohesion array
    fill!(afX,zero(T))
    fill!(afY,zero(T))

    for tp=1:length(pList), tw=1:length(wList)
        NearestPoint!(pointOnWall,pList[tp],wList[tw])

        if isCloseEnough(pointOnWall,pList[tp],k,rList[tp])
            generateQuadNodes!(xquad,yquad,pointOnWall,pList[tp],wList[tw],k,rList[tp])
            Fx,Fy = wallTrapQuad(pList[tp],wList[tw],xquad,yquad,ϵ,rList[tp],k,pointOnWall)
            afX[tp] += Fx
            afY[tp] += Fy
        end 
    end

    nothing
end