module eFEMpart

using Reexport

@reexport using eFEM
@reexport using StokesParticles
@reexport using figtree_jll

#include("fem_sp_coupling.jl")

end # module 
