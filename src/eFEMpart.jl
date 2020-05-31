module eFEMpart

using Reexport

@reexport using eFEM
@reexport using StokesParticles

export update_particles_noFEMcl!

export scratch_data

include("struct_defs.jl")
include("barycentric_interpolation.jl")
include("seepage_force.jl")
include("update_particles.jl")

#include("fem_sp_coupling.jl")

end # module 
