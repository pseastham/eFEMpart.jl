# File containing functions to update particle positions based on different forces

# computes seepage velocities at particle positions
include("SeepageForce.jl")   # itself calls CellLists.jl, isInside.jl, BarycentricInterpolation.jl
include("CohesionForce.jl")  # itself calls CellLists.jl
include("AdhesionForce.jl")  # itself calls CellLists.jl, SinkholeTypes.jl, CohesionForce.jl, isInside.jl

# Uses:
#   1. Graviational Force
#   2. Seepage Force (fluid velocity)
#   3. particle-particle Cohesion Force (Lennard Jones Potential)
#   4. particle-particle Adhesion Force (Lennard Jones Potential)
function updateParticle_all!(mesh,pList::Vector{Point2D},rList::Vector{T},wList::Vector{W},
                                u::Vector{T},v::Vector{T},paramArr::Vector{Any},Δt::T) where {T<:Real, W<:AbstractWall}
    pUarr,pVarr = computeParticleVelocity_all(mesh,pList,rList,wList,u,v,paramArr)

    # =========================
    # Update particle position
    # =========================
    for ti=1:length(rList)
        pList[ti].x += pUarr[ti]*Δt
        pList[ti].y += pVarr[ti]*Δt
    end

    nothing
end

"""
    ComputeParticleVelocity_all!

Computes particle velocity for all particles. Goal is to have zero memory allocation and be completely type-stable

INPUT:
    mesh:      finite element mesh
    pList:     list of particle positions (in Point2D datatype)
    rList:     list of particle radii
    wList:     list of all walls for particles to interact against
    u,v:       2D velocity components 
    paramArr:  array of scalar parameters -- might switch out for Parameters.jl datatypes

OUTPUT:
    // nothing (updates BLANK in-place)
    pUarr, pVarr: 2D velocity arrays for particles
"""
function computeParticleVelocity_all(mesh,pList::Vector{Point2D},rList::Vector{T},wList::Vector{W},
                                        u::Vector{T},v::Vector{T},paramArr::Vector{Any}) where {T<:Real, W<:AbstractWall}
    # G:      gravitation parameter
    # rc:     cutoff ratio
    # ϵ:      strength of LJ force   
    # k:      multiple of equilibrium below which we should consider adhesion force calculations (i.e. k=2, k=5, etc)
    # Nquad:  number of quadrature nodes for particle-wall interaction (adhesion)
    G,rc,ϵ,k,Nquad = paramArr

    Nparticles = length(pList)

    # -----------------
    # USED FOR TESTING
    # -----------------
    particleVolume  = 1.0
    particleDensity = 1.0

    # ----------------------------------------------------------------
    # INITIALIZATIONS -- EVENTUALLY WILL BE MOVED OUTSIDE OF FUNCTION
    # initialize arrays to don't need to re-allocate memory unnecessarily during velocity interpolation
    # ----------------------------------------------------------------
    # arrays for barycentric interpolation
    polygon      = [Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0)]
    extremePoint = Point2D(100_000.0,0.0)
    w   = zeros(4)           
    uEl = zeros(4)
    vEl = zeros(4)
    a = zeros(3)
    b = zeros(3)
    c = zeros(3)
    d = zeros(3)
    pointOnWall = Point2D(0.0,0.0)
    xquad = zeros(Float64,Nquad)
    yquad = zeros(Float64,Nquad)

    # seepage velocities -- will be updated in-place
    uSeepage = zeros(Float64,Nparticles)
    vSeepage = zeros(Float64,Nparticles)

    # cohesion force -- will be updated in-place
    cfX = zeros(Float64,Nparticles)
    cfY = zeros(Float64,Nparticles)

    # adhesion force -- will be updated in-place
    afX = zeros(Float64,Nparticles)
    afY = zeros(Float64,Nparticles)

    # GENERATE PARTICLE CELL LIST --
    # in future might want to create 2 particle cell lists,
    # one for seepage calculation and one for cohesion calculation
    # where clL changes between the two for increased efficiency
    clL = 0.1
    clTotalBounds = [-1.0,1.0,-1.0,1.0]         # WILL LIKELY NEED TO CHANGE!
    particleCL = generateCellList(pList,clTotalBounds,clL)

    # generate cell list - FEM map 
    femCLmap = FEMgenerateMap(mesh,clTotalBounds,clL)

    # 1. compute gravitational force
    gfX = zeros(Nparticles) 
    gfY = -G*ones(Nparticles)

    # 2. interpolate seepage velocity of fluid at position of particles
    # REQUIRES CELL LIST for particles and meshCLmap for FEM <-> cell list indexing of elements
    # compute cohesion force on all particles using cell list
    BarycentricVelocityInterp_CL!(uSeepage,vSeepage,mesh,u,v,pList,polygon,particleCL,
                                    femCLmap,extremePoint,w,uEl,vEl,a,b,c,d)

    # 3. compute cohesion forces
    computeCohesion_CL!(cfX,cfY,pList,rList,rc,ϵ,particleCL)

    # 4. compute adhesion forces
    AdhesionForce!(afX,afY,pList,rList,wList,k,rc,ϵ,pointOnWall,xquad,yquad)

    # 5. use stokes force balance to compute particle velocities
    #κ = permeability
    κ = 1.0
    #n = porosity
    n = 1.0
    #gammas = specific weight
    gammas = 1.0
    pUarr = κ*(gfX + cfX + afX + n*gammas*uSeepage/κ)/(n*gammas)
    pVarr = κ*(gfY + cfY + afY + n*gammas*vSeepage/κ)/(n*gammas)

    return pUarr,pVarr
end