using SafeTestsets

@time begin
    using eFEMpart
#    @time @safetestset "Cell List Tests" begin include("cl_test.jl") end
end
