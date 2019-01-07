# This module loads functions and types necessary for driver files
#    used in Sinkhole problems

module SinkholeModule

# loads all finite element stuff
using eFEM

# SinkholeTypes.jl
export AbstractWall,            # Abstract type that encapsulates all concrete wall types
       LineWall,                # mutable type for describing line walls
       CircleWall,              # mutable type for describing circle walls
       ArcWall,                 # mutable type for describing circle arc walls
       Particle,
       Point,                   # mutable type for points
       ParticleList             # Initializes list of particles

# CellLists.jl
export CellGenerate             # Generates cell lists from either a list
                                #   of particles or a mesh

# Porosity.jl
export SandToPorosityGaussian

# ParticleForces.jl
export LennardJonesForce,
       LennardJonesPotential,
       LennardJonesPotentialMagnitude

# WallRepulsion.jl 
export NearestPoint,
       isCloseEnough,
       GenerateQuadNodes,
       isInLine,
       WallTrapQuad

# UpdateParticles.jl
export UpdateParticle!,
       UpdateParticle_Gravity!,
       UpdateParticle_LJWalls!,
       UpdateParticle_Cohesion!,
       UpdateParticle_Seepage!

# include function files
include("SinkholeTypes.jl")
include("Walls.jl")
include("CellLists.jl")
include("Porosity.jl")
include("ParticleForces.jl")
include("WallRepulsion.jl")
include("UpdateParticles.jl")

end # SinkholeModule
