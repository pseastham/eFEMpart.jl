# File to be loaded into eFEMpart

# Functions for computing contact forces of 
# particles and walls via LJ potentials
# and a quadrature

### Compute nearest point to particle on LineWall, CircleWall, and ArcWall
function NearestPoint(particle::Particle,wall::LineWall)
    px=particle.xpos; py=particle.ypos

    Ax=wall.nodes[1].x; Ay=wall.nodes[1].y
    Bx=wall.nodes[2].x; By=wall.nodes[2].y

    bx=px-Ax; by=py-Ay
    ax=Bx-Ax; ay=By-Ay

    ℓ = sqrt(ax^2+ay^2)

    dotprod = ax*bx + ay*by
    
    sx = dotprod*ax/ℓ^2 + Ax
    sy = dotprod*ay/ℓ^2 + Ay

    return [sx;sy]
end
function NearestPoint(particle::Particle,wall::CircleWall)
    px=particle.xpos; py=particle.ypos
    cx=wall.center.x; cy=wall.center.y
    r = wall.radius
    θ = atan2(py-cy,px-cx)

    sx = cx + r*cos(θ)
    sy = cy + r*sin(θ)

    return [sx;sy]
end
function NearestPoint(particle::Particle,wall::ArcWall)
    px= particle.xpos; py=particle.ypos
    cx= wall.nodes[2].x; cy=wall.nodes[2].y
    r = sqrt((cx-wall.nodes[1].x)^2 + (cy-wall.nodes[1].y)^2)
    θ = atan2(py-cy,px-cx)

    sx = cx + r*cos(θ)
    sy = cy + r*sin(θ)

    return [sx;sy]
end

function isCloseEnough(p::Particle,s::Vector{Float64},k::Float64,rm::Float64)
    px = p.xpos; py = p.ypos
    ℓ = sqrt((px-s[1])^2 + (py-s[2])^2)

    return  (ℓ < k*rm ? true : false)
end

function GenerateQuadNodes(p::Particle,wall::LineWall,s::Vector{Float64},k::Float64,rm::Float64,N::Int)
    Ax=wall.nodes[1].x; Ay=wall.nodes[1].y
    tx=wall.t[1]; ty=wall.t[2]

    projLen = sqrt((p.xpos - s[1])^2 + (p.ypos-s[2])^2)
    L = 2*rm

    Δs = 2*L/(N-1)
    xquad = [s[1] - L*tx + (i-1)*Δs*tx for i=1:N]
    yquad = [s[2] - L*ty + (i-1)*Δs*ty for i=1:N]

    return xquad, yquad
end
function GenerateQuadNodes(p::Particle,wall::CircleWall,s::Vector{Float64},k::Float64,rm::Float64,N::Int)
    px=p.xpos; py=p.ypos
    cx=wall.center.x; cy=wall.center.y
    r = wall.radius
    θcenter = atan2(py-cy,px-cx)
    θmin = θcenter - pi/2
    Δθ = pi/(N-1)

    xquad = [cx + r*cos(θmin+(i-1)*Δθ) for i=1:N]
    yquad = [cy + r*sin(θmin+(i-1)*Δθ) for i=1:N]

    return xquad, yquad
end
function GenerateQuadNodes(p::Particle,wall::ArcWall,s::Vector{Float64},k::Float64,rm::Float64,N::Int)
    px=p.xpos; py=p.ypos
    cx= wall.nodes[2].x; cy=wall.nodes[2].y
    r = sqrt((cx-wall.nodes[1].x)^2 + (cy-wall.nodes[1].y)^2)
    θcenter = atan2(py-cy,px-cx)
    θmin = θcenter - pi/2
    Δθ = 2*pi/(N-1)

    xquad = [cx + r*cos(θmin+(i-1)*Δθ) for i=1:N]
    yquad = [cy + r*sin(θmin+(i-1)*Δθ) for i=1:N]

    return xquad, yquad
end

# function to determine whether quadrature node is within line
function isInLine(wall::LineWall,sx,sy)
    Ax=wall.nodes[1].x; Ay=wall.nodes[1].y
    Bx=wall.nodes[2].x; By=wall.nodes[2].y
    A=[Ax,Ay]; B=[Bx,By]; s=[sx,sy]
    return onSegment(A,s,B)
end
function isInLine(wall::CircleWall,sx,sy)
    cx = wall.center.x; cy=wall.center.y
    val = abs(wall.radius - sqrt((cx-sx)^2 + (cy-sy)^2))
    TOL = 1e-12
    if val < TOL
        return true
    else
        return false
    end
end
function isInLine(wall::ArcWall,sx,sy)
    cx = wall.nodes[2].x; cy=wall.nodes[2].y
    radius = sqrt((cx-wall.nodes[1].x)^2 + (cy-wall.nodes[1].y)^2)
    val = abs(radius - sqrt((cx-sx)^2 + (cy-sy)^2))

    θ1 = atan2(wall.nodes[1].y-cy,wall.nodes[1].x-cx)
    θ2 = atan2(wall.nodes[3].y-cy,wall.nodes[3].x-cx)
    θs = atan2(sy-cy,sx-cx)

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
# ϵ = strength of lennard jones potential
# rm = equib distance for lennard jones porential
function WallTrapQuad(p::Particle,wall::LineWall,
                      xquad::Vector{Float64},yquad::Vector{Float64},
                      ϵ::Float64,rm::Float64)
    N = length(xquad)
    Δs = sqrt((wall.nodes[1].x - wall.nodes[2].x)^2 + (wall.nodes[1].y - wall.nodes[2].y)^2)/(N-1)
    
    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=2:N
        if isInLine(wall,xquad[ti],yquad[ti])
            # compute vector from x,y to p 
            qxKm1 = p.xpos-xquad[ti-1]; qyKm1 = p.ypos-yquad[ti-1]
            qxK   = p.xpos-xquad[ti];   qyK   = p.ypos-yquad[ti]

            # compute strength of LJ force
            FxKm1,FyKm1 = LennardJonesForce(ϵ,rm,qxKm1,qyKm1)
            FxK,FyK     = LennardJonesForce(ϵ,rm,qxK,qyK)

            quadSumX += 0.5*(FxKm1 + FxK)*Δs
            quadSumY += 0.5*(FyKm1 + FyK)*Δs
        end
    end
    return quadSumX,quadSumY
end
function WallTrapQuad(p::Particle,wall::CircleWall,
                    xquad::Vector{Float64},yquad::Vector{Float64},
                    ϵ::Float64,rm::Float64)
    cx = wall.center.x; cy = wall.center.y
    θ1 = atan2(yquad[1]-cy,xquad[1]-cx)
    θ2 = atan2(yquad[2]-cy,xquad[2]-cx)
    Δs = abs(wall.radius*(θ2-θ1))
    N = length(xquad)
    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=1:N
        if isInLine(wall,xquad[ti],yquad[ti])
            # compute vector from x,y to p 
            qx=p.xpos - xquad[ti]; qy=p.ypos - yquad[ti]

            # compute strength of LJ force
            Fx,Fy = LennardJonesForce(ϵ,rm,qx,qy)

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
function WallTrapQuad(p::Particle,wall::ArcWall,
                    xquad::Vector{Float64},yquad::Vector{Float64},
                    ϵ::Float64,rm::Float64)
    cx = wall.nodes[2].x; cy = wall.nodes[2].y
    θ1 = atan2(yquad[1]-cy,xquad[1]-cx)
    θ2 = atan2(yquad[2]-cy,xquad[2]-cx)
    radius = sqrt((cx - wall.nodes[1].x)^2 + (cy - wall.nodes[1].y)^2)

    if θ1 > θ2
        θ2 += 2pi
    end

    Δs = abs(radius*(θ2-θ1))
    N = length(xquad)
    quadSumX = 0.0
    quadSumY = 0.0
    # for each particle
    for ti=1:N
        if isInLine(wall,xquad[ti],yquad[ti])
            # compute vector from x,y to p 
            qx=p.xpos - xquad[ti]; qy=p.ypos - yquad[ti]

            # compute strength of LJ force
            Fx,Fy = LennardJonesForce(ϵ,rm,qx,qy)

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