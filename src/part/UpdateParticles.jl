# Code to be loaded into eFEMpart

# File containing functions to update particle positions based on different forces

# Uses:
#   1. Graviational Force
#   2. Seepage Force
#   3. Cohesion Force (Lennard Jones Potential)
#
# OBSOLETE!!!!
function UpdateParticle!(particle,j::Int,dt::Float64,
                         pList,wList,cl,rm::Float64,FEMmesh,solution,square)
    # ==============================
    # 1. Compute Graviational Force
    # ==============================
    g = 1.0
    FxG,FyG = GravitationForce(g)

    # =========================
    # 2. compute Seepage velocity
    # =========================
    if isInsideRect(square,[particle.xpos,particle.ypos])
        # interpolate velocity at particle from FEM solution
        fluidVel = FEMinterpolateVelocity(particle,FEMmesh,solution)
    else
        fluidVel = [0.0,0.0]
    end

    # ===========================================
    # 3. Compute Cohesion Force using Cell Lists
    # ===========================================
    FxC,FyC = CohesionForceCL(particle,j,pList,rm,cl)

    # ========================
    # Sum all forces together
    # ========================
    # need to interpolate the porosity/permeability
    #κ = permeability
    κ = 1.0
    #n = porosity
    n = 1.0
    #gammas = specific weight
    gammas = 1.0
    partVelx = κ*(FxG + FxC + n*gammas*fluidVel[1]/κ)/(n*gammas)
    partVely = κ*(FyG + FyC + n*gammas*fluidVel[2]/κ)/(n*gammas)

    # ========================================================
    # checks for boundary, and if so, bounces particle off it
    # ========================================================
    bounce,wIndex = checkInsideWall(particle,wList)
    if bounce
    partVelx,partVely = WallBounce(wList[wIndex],[partVelx,partVely])
    end

    # =========================
    # Update particle position
    # =========================
    particle.xpos += partVelx*dt
    particle.ypos += partVely*dt
end

function UpdateParticle_Gravity!(particle::Particle,Δt,G,κ,n,γs)
    # compute force of gravity
    FxG,FyG = GravitationForce(G)

    # Stokes Force Balance to compute velocity
    particleVelocityX = κ*FxG/(n*γs)
    particleVelocityY = κ*FyG/(n*γs)

    # Update particle position
    particle.xpos += particleVelocityX*Δt
    particle.ypos += particleVelocityY*Δt
end

function UpdateParticle_LJWalls!(particle::Particle,wList::Vector{T},Δt,κ,n,γs,
                                k::Float64,rm::Float64,ϵ::Float64,Nint::Int) where T<:AbstractWall
    Nwalls = length(wList)

    # sum up total forces in x and y directions
    FxW = 0.0; FyW = 0.0
    for i=1:Nwalls
        Fxtemp,Fytemp = LJWallForce(particle,wList[i],k,rm,ϵ,Nint)
        FxW += Fxtemp; FyW += Fytemp
    end

    # Stokes Force Balance to compute velocity
    particleVelocityX = κ*FxW/(n*γs)
    particleVelocityY = κ*FyW/(n*γs)

    # Update particle position
    particle.xpos += particleVelocityX*Δt
    particle.ypos += particleVelocityY*Δt
end

function UpdateParticle_Cohesion!(particle::Particle,j::Int,pList,rm,ϵ,cl::Vector{Cell},Δt,κ,n,γs)
    # Compute Cohesion Force using Cell Lists
    FxC,FyC = CohesionForceCL(particle,j,pList,ϵ,rm,cl)

    # Stokes Force Balance to compute velocity
    particleVelocityX = κ*FxC/(n*γs)
    particleVelocityY = κ*FyC/(n*γs)

    # Update particle position
    particle.xpos += particleVelocityX*Δt
    particle.ypos += particleVelocityY*Δt
end

function UpdateParticle_Seepage!(particle::Particle,FEMmesh,DarcySolution,square,
                                 Δt,κ,n,γs)
    # Compute Seepage Force using Cell Lists
    if isInsideRect(square,[particle.xpos,particle.ypos])
        # interpolate velocity at particle from FEM solution
        FxS,FyS = FEMinterpolateVelocity(particle,FEMmesh,DarcySolution)
    else
        FxS,FyS = [0.0,0.0]
    end

    # Stokes Force Balance to compute velocity
    particleVelocityX = κ*FxS/(n*γs)
    particleVelocityY = κ*FyS/(n*γs)

    # Update particle position
    particle.xpos += particleVelocityX*Δt
    particle.ypos += particleVelocityY*Δt
end

function UpdateParticle_All!(particle::Particle,j::Int,FEMmesh,pList,wList::Vector{T},Δt,G,κ,n::Float64,γs,
                            k,rm,ϵp,ϵw,Nint::Int,cl::Vector{Cell},
                            DarcySolution,square) where T<:AbstractWall

end