using BenchmarkTools

include("../../src/part/BarycentricInterpolation.jl")

function computeBaryWeights_TEST()
    x = [0.0, 1.0, 0.8, 0.0]
    y = [0.0, 0.0, 0.8, 0.4]
    p = [0.8, 0.2]

    w = computeBaryWeights(x, y, p)

    @btime computeBaryWeights($x,$y,$p)

    z = [0.0, 0.5, 1.0, 1.5]
    zi = dot(z, w)
end
