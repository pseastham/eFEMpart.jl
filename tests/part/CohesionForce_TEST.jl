# computes cohesion forces

using BenchmarkTools

include("../../src/part/CellLists.jl")   # for Point2D

function LennardJonesForceMagnitude_TEST()
    ϵ  = 1.0
    d  = 1.0
    r  = 1.0 
    rc = 2.0

    #@code_warntype LennardJonesForceMagnitude(ϵ,d,r,rc)
    @benchmark LennardJonesForceMagnitude($ϵ,$d,$r,$rc)

    #nothing
end

function ForceCalculation_TEST()
    ϵ  = 1.0
    d  = 1.0
    r  = 1.0 
    rc = 1.0

    Δx = 2.0
    Δy = 1.5

    #@code_warntype ForceCalculation(ϵ,d,rc,Δx,Δy)
    @benchmark ForceCalculation($ϵ,$d,$rc,$Δx,$Δy)

    #nothing
end

# N: number of particles
function computeCohesion_TEST(N=100)
    nodeList = [Point2D(2*rand()-1.0,2*rand()-1.0) for i=1:N]
    radiusList = ones(N)

    #@code_warntype computeCohesion(nodeList,radiusList)
    @btime computeCohesion_SMARTER($nodeList,$radiusList)

    nothing
end

function computeCohesion_CL!_TEST(N=100,L=0.1)
    nodeList = [Point2D(2*rand()-1.0,2*rand()-1.0) for i=1:N]
    radiusList = 0.5*L*ones(N)

    # construct cell list
    totalBounds = [-1.0-L,1.0+L,-1.0-L,1.0+L]
    cl = generateCellList(nodeList,totalBounds,L)

    # initialize vector for force assignment
    cfX = zeros(Float64,N)
    cfY = zeros(Float64,N)

    #@code_warntype computeCohesion_CL!(nodeList,radiusList,cl,cfX,cfY)
    @btime computeCohesion_CL!($nodeList,$radiusList,$cl,$cfX,$cfY)

    nothing
end

# computes L (cell list side length) ideal for lennard-jones force calculations
# from particle-particle interactions (cohesion)
function computeLJ_L(rc::T,rList::Vector{T}) where T<:Real
    d = 2*maximum(rList)
    return 1.1*rc*d/sqrt(2)
end
function computeLJ_L(rc::T,radius::T) where T<:Real
    d = 2*radius
    return 1.1*rc*d/sqrt(2)
end

# ---------------
# ACCURACY TESTS
# ---------------

function force_plot_TEST()
    ϵ  = 1.0         # strength
    rc = 1.5         # ratio to cut-off

    N = 30
    L = 0.1
    #nodeList = [Point2D(-L+L/3,0.0),Point2D(L-L/3,0.0)]

    nodeList = [Point2D(2*rand()-1.0,2*rand()-1.0) for i=1:N]
    radiusList = 0.5*L*ones(N)

    # initalize force arrays
    cfX = zeros(N); cfY = zeros(N)

    # direct
    computeCohesion!(nodeList,radiusList,rc,ϵ,cfX,cfY)

    for i=1:N
        println(i,"\t",sqrt(cfX[i]^2 + cfY[i]^2))
    end


    return nodeList,radiusList,cfX,cfY,rc
end

function twoParticleComparison_TEST(N)
    ϵ  = 1.0         # strength
    rc = 1.5         # ratio to cut-off

    nodeList = [Point2D(2*rand()-1.0,2*rand()-1.0) for i=1:N]
    radiusList = 0.01*ones(N)

    L = computeLJ_L(rc,radiusList)
    #println("L: ",L)

    # check that minimum cell list length is satisfied
    if (L < rc*2*radiusList[1]/sqrt(2))
        println("L:   ",L)
        println("min: ",rc*2*radiusList[1]/sqrt(2))
        throw(error("cell list length must be larger!"))
    end


    # construct cell list
    totalBounds = [-1.0-L,1.0+L,-1.0-L,1.0+L]
    cl = generateCellList(nodeList,totalBounds,L)

    # initalize force arrays
    cfXdirect = zeros(N); cfYdirect = zeros(N)
    cfXcl     = zeros(N); cfYcl     = zeros(N)

    # direct
    println("N: ",N)
    print("direct:    ")
    #@code_warntype computeCohesion!(nodeList,radiusList,rc,ϵ,cfXdirect,cfYdirect)
    @btime computeCohesion!($nodeList,$radiusList,$rc,$ϵ,$cfXdirect,$cfYdirect)

    # cell list
    print("cell list: ")
    #@code_warntype computeCohesion_CL!(nodeList,radiusList,rc,ϵ,cl,cfXcl,cfYcl)
    @btime computeCohesion_CL!($nodeList,$radiusList,$rc,$ϵ,$cl,$cfXcl,$cfYcl)
    println()

    if false
        xdiff = zeros(N)
        ydiff = zeros(N)
        for i=1:N
            if cfXdirect[i] != 0.0
                xdiff[i] = abs((cfXdirect[i] - cfXcl[i])/cfXdirect[i])
            else
                xdiff[i] = abs(cfXdirect[i] - cfXcl[i])
            end

            if cfYdirect[i] != 0.0
                ydiff[i] = abs((cfYdirect[i] - cfYcl[i])/cfYdirect[i])
            else
                ydiff[i] = abs(cfYdirect[i] - cfYcl[i])
            end
        end

        if any(xdiff .> 1e-12) || any(ydiff .> 1e-12)
            println("force vectors are different!")
            println("max difference: ",maximum([xdiff;ydiff]))
        else
            println("force vectors match!")
        end

        println()

        println("X: ")
        for i=1:N
            println(i,"\t cfXd:  ",cfXdirect[i])
            println("\t cfXcl: ",cfXcl[i])
        end
        println()
        println("Y: ")
        for i=1:N
            println(i,"\t cfYd:  ",cfYdirect[i])
            println("\t cfYcl: ",cfYcl[i])
        end
    end

    nothing
end

function main()
    for N in [100,200,400,800,1600,3200,6400,12800]
        twoParticleComparison_TEST(N);
    end

    nothing
end