# File to be loaded into eFEMpart

# This file contains functions used to compute forces on particles

# updates position of particle by computing velocity,
# uses cell list to minimize cohesion force computation
# and Stokes Force Balance (not F=ma)

# ==============================================================================
# FORCES
# ==============================================================================

# returns vector representing force of gravity in -y direction
function GravitationForce(g::Float64)
    return 0.0,-g
end

function BuoyancyForce(pressureGradient)
  porosity = 1.0
  SWS = 1.0             # specific weight of sand
  n = 1.0

  Fx = -(1-n)*pressureGradient[1]
  Fy = -(1-n)*pressureGradient[2]

  return Fx,Fy
end

# returns seepage force
function SeepageForce(particle,velocity)
  porosity = 1.0
  SWS = 1.0             # specific weight of sand
  permeability = 1.0

  Fx = porosity*SWS*velocity[1]/permeability
  Fy = porosity*SWS*velocity[2]/permeability

  return Fx,Fy
end

function FEMinterpolateVelocity(particle,FEMmesh,solution)
    u = pointTransform(FEMmesh,solution.u,[particle.xpos, particle.ypos])
    v = pointTransform(FEMmesh,solution.v,[particle.xpos, particle.ypos])

    return u,v
end

function FEMinterpolateScalar(particle,FEMmesh,solution)
    return pointTransform(FEMmesh,solution.u,[particle.xpos, particle.ypos])
end

# returns the lennard-jones potential V(r)
function LennardJonesPotential(ϵ,rm,r)
    ratio = rm/r
    f1 = ratio^12
    f2 = ratio^6

    return ϵ*(f1 -2*f2)
end

# Returns -d/dr V(r) where V(r) is the lennard-jones potential
function LennardJonesPotentialMagnitude(ϵ::T,rm::T,r::T) where T<:Real
    if r < 0.96*rm
        r = 0.96*rm
        return 12*ϵ*rm^6*(rm-r)*(rm+r)*(rm^2 + r^2 + rm*r)*(rm^2 + r^2 - rm*r)/r^13
    elseif r > 1.9*rm 
        return 0.0
    else
        return 12*ϵ*rm^6*(rm-r)*(rm+r)*(rm^2 + r^2 + rm*r)*(rm^2 + r^2 - rm*r)/r^13
    end
end

function LennardJonesForce(ϵ::T,rm::T,Δx::T,Δy::T) where T<:Real
    len = sqrt(Δx^2 + Δy^2)
    tx = Δx/len
    ty = Δy/len

    LJmag = LennardJonesPotentialMagnitude(ϵ,rm,len)

    Fx = tx*LJmag
    Fy = ty*LJmag

    return Fx,Fy
end

function CohesionForceCL(particle,j,pList,ϵ,rm,cl)
    Npl = length(pList)
    Ncl = length(cl)
    pcc = 0             # cell index containing particle

    # find cell containing particle
    for i=1:Ncl
        if j in cl[i].ParticleList
            pcc = i
            break
        end
    end

    if pcc !=0
        # find neighbors of particle cell
        nCells = vcat(pcc,cl[pcc].NeighborList)

        # find particles in all relevant cells
        pclose = Array{Int,1}(undef,0)
        for i=1:length(nCells)
            ncl = length(cl[nCells[i]].ParticleList)
            if ncl > 0
                for tj=1:ncl
                    push!(pclose,cl[nCells[i]].ParticleList[tj])
                end
            end
        end

        x1 = particle.xpos
        y1 = particle.ypos

        Fx = 0.0
        Fy = 0.0
        for tj in pclose
            Δx = x1 - pList[tj].xpos
            Δy = y1 - pList[tj].ypos

            # avoid self-particle
            if abs(Δx+Δy)>1e-12
                dFx,dFy = LennardJonesForce(ϵ,rm,Δx,Δy)
                Fx += dFx
                Fy += dFy
            end
        end
    else
        Fx = 0.0
        Fy = 0.0
    end

    return Fx,Fy
end

function LJWallForce(particle::Particle,wall::T,k::Float64,rm::Float64,ϵ::Float64,Nint::Int) where T<:AbstractWall
    s = NearestPoint(particle,wall)
    close = isCloseEnough(particle,s,k,rm)

    Fx = 0.0; Fy = 0.0
    if close
        xquad,yquad = GenerateQuadNodes(particle,wall,s,k,rm,Nint)

        # compute force
        Fx,Fy = WallTrapQuad(particle,wall,xquad,yquad,ϵ,rm)
    end

    return Fx,Fy
end
