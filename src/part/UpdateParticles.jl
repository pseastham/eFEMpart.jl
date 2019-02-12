# Code to be loaded into eFEMpart

# File containing functions to update particle positions based on different forces

# Uses:
#   1. Graviational Force
#   2. Seepage Force (velocity)
#   3. Cohesion Force (Lennard Jones Potential)
#   4. Wall Force (Lennard Jones Potential)
function UpdateParticle_all!(particle,j::Int,pList,wList,cl,rm::Float64,
                             Δt,k,ϵp,ϵw,Nint,square,FEMmesh,solution,
                             G,κ,n,gammas)
    # ==============================
    # 1. Compute Graviational Force
    # ==============================
    FxG,FyG = GravitationForce(G)

    # ============================
    # 2. compute Seepage velocity
    # ============================
    if isInsideRect(square,[particle.xpos,particle.ypos])
        # interpolate velocity at particle from FEM solution
        fluidVel = FEMinterpolateVelocity(particle,FEMmesh,solution)
    else
        fluidVel = [0.0,0.0]
    end
    # FOR TESTING!!!
    #fluidVel = [0.0,0.0]

    # ===========================================
    # 3. Compute Cohesion Force using Cell Lists
    # ===========================================
    #FxC,FyC = CohesionForceCL(particle,j,pList,ϵp,rm,cl)
    FxC,FyC = 0.0,0.0

    # ===========================================
    # 4. Compute Wall Force
    # ===========================================
    Nwalls = length(wList)

    # sum up total forces in x and y directions
    FxW = 0.0; FyW = 0.0
    for i=1:Nwalls
        Fxtemp,Fytemp = LJWallForce(particle,wList[i],k,rm,ϵw,Nint)
        FxW += Fxtemp; FyW += Fytemp
    end

    # ========================
    # Sum all forces together
    # ========================
    # need to interpolate the porosity/permeability
    # κ = permeability
    # n = porosity
    # gammas = specific weight
    partVelx = κ*(FxG + FxC + FxW + n*gammas*fluidVel[1]/κ)/(n*gammas)
    partVely = κ*(FyG + FyC + FyW + n*gammas*fluidVel[2]/κ)/(n*gammas)

    # =========================
    # Update particle position
    # =========================
    particle.xpos += partVelx*Δt
    particle.ypos += partVely*Δt
end

# This one is to be used without velocity
function UpdateParticle_novelocity!(particle,j::Int,pList,wList,cl,rm::Float64,
                                    Δt,k,ϵp,ϵw,Nint,G)
    # ==============================
    # 1. Compute Graviational Force
    # ==============================
    FxG,FyG = GravitationForce(G)

    # =========================
    # 2. compute Seepage velocity -- no velocity!
    # =========================
    fluidVel = [0.0,0.0]

    # ===========================================
    # 3. Compute Cohesion Force using Cell Lists
    # ===========================================
    FxC,FyC = CohesionForceCL(particle,j,pList,ϵp,rm,cl)


    # ===========================================
    # 4. Compute Wall Force
    # ===========================================
    Nwalls = length(wList)

    # sum up total forces in x and y directions
    FxW = 0.0; FyW = 0.0
    for i=1:Nwalls
        Fxtemp,Fytemp = LJWallForce(particle,wList[i],k,rm,ϵw,Nint)
        FxW += Fxtemp; FyW += Fytemp
    end

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
    partVelx = κ*(FxG + FxC + + FxW + n*gammas*fluidVel[1]/κ)/(n*gammas)
    partVely = κ*(FyG + FyC + + FyW + n*gammas*fluidVel[2]/κ)/(n*gammas)

    # =========================
    # Update particle position
    # =========================
    particle.xpos += partVelx*Δt
    particle.ypos += partVely*Δt
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