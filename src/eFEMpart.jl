module eFEMpart

using Reexport

@reexport using eFEM
@reexport using StokesParticles

include("fem_sp_coupling.jl")

end # module 
