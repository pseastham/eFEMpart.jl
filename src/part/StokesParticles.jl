module StokesParticles

export AbstractWall, NearestPoint, generateQuadNodes!,isInLine

# CellLists.jl
export femGenerateMap,
       generateCellList,
       updateCellList!

# fgt.jl
export interpFGT!

# ParticleTypes.jl
export Point2D,
       LineWall,
       CircleWall,
       ArcWall

# UpdateParticles.jl
export updateParticle_all!,
       updateParticle_all_nofluid!

# for testing...
export BarycentricVelocityInterp_CL!,
       computeCohesion!,
       AdhesionForce!,
       computeCohesion_CL!,
       LennardJonesForceMagnitude

include("part/ParticleTypes.jl")
include("part/CellLists.jl")
include("part/fgt.jl")
include("part/UpdateParticles.jl")

include("part/AdhesionForce.jl")

end # module
