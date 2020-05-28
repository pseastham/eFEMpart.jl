using SafeTestsets

@time begin
    @time @safetestset "Cell List Tests" begin include("cl_test.jl") end
end
