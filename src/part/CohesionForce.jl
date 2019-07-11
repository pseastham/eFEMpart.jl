# computes cohesion forces

#using eFEMpart
using BenchmarkTools

# for Point2D
include("CellLists.jl")

# Returns -d/dr V(r) where V(r) is the lennard-jones potential
# ϵ: strength
# d: equilibrium distance
# r: distance between particles
# rc: cutoff distance as multiple of d 
function LennardJonesForceMagnitude(ϵ::T,d::T,r::T,rc::T) where T<:Real
    if r > rc*d
        return zero(T)
    else
        d2 = d*d
        r2 = r*r
        dr = d*r
        return -12*ϵ*d^6*(d-r)*(d+r)*(d2 + r2 + dr)*(d2 + r2 - dr)/r^13
    end
end

function ForceCalculation(ϵ::T,d::T,rc::T,Δx::T,Δy::T) where T<:Real
    len = sqrt(Δx*Δx + Δy*Δy)
    tx=Δx/len; ty=Δy/len

    # lennard-jones potential
    LJmag = LennardJonesForceMagnitude(ϵ,d,len,rc)

    Fx = tx*LJmag
    Fy = ty*LJmag

    return Fx,Fy
end 

function computeCohesion!(cfX::Vector{T},cfY::Vector{T},nodeList::Vector{Point2D},radiusList::Vector{T},rc::T,ϵ::T) where T<:Real
    Nparticles = length(nodeList)

    # zero-out cohesion
    fill!(cfX,zero(T))
    fill!(cfY,zero(T))

    for ti=1:Nparticles, tj=(ti+1):Nparticles
        Δx = nodeList[tj].x - nodeList[ti].x
        Δy = nodeList[tj].y - nodeList[ti].y

        d = radiusList[ti] + radiusList[tj]

        fx,fy = ForceCalculation(ϵ,d,rc,Δx,Δy)
        cfX[ti] += fx
        cfY[ti] += fy
        cfX[tj] -= fx
        cfY[tj] -= fy
    end

    nothing
end

# nodelist: list of points to compute force between
# radiuslist: list of particle radii            (currently unused)
# cfX: cohesion force in x-direction
# cfY: cohesion force in y-direction
function computeCohesion_CL!(cfX::Vector{T},cfY::Vector{T},nodeList::Vector{Point2D},radiusList::Vector{T},
                                rc::T,ϵ::T,cl::CellList) where T<:Real
    # zero-out cohesion
    fill!(cfX,zero(T))
    fill!(cfY,zero(T))

    # go over cell lists
    for cellInd = 1:length(cl.cells)
        # compute cell-cell forces
        numNodes = length(cl.cells[cellInd].nodeList)
        for ti=1:numNodes, tj=(ti+1):numNodes
            nodeI = cl.cells[cellInd].nodeList[ti]
            nodeJ = cl.cells[cellInd].nodeList[tj]

            Δx = nodeList[nodeJ].x - nodeList[nodeI].x
            Δy = nodeList[nodeJ].y - nodeList[nodeI].y
    
            d = radiusList[nodeI] + radiusList[nodeJ]
    
            fx,fy = ForceCalculation(ϵ,d,rc,Δx,Δy)
            cfX[nodeI] += fx; cfY[nodeI] += fy
            cfX[nodeJ] -= fx; cfY[nodeJ] -= fy
        end

        # compute cell-neighbor forces
        for ncInd in cl.cells[cellInd].neighborList
            numNeighborNodes = length(cl.cells[ncInd].nodeList)
            for ti=1:numNodes, tj=1:numNeighborNodes
                nodeI = cl.cells[cellInd].nodeList[ti]
                nodeJ = cl.cells[ncInd].nodeList[tj]

                Δx = nodeList[nodeJ].x - nodeList[nodeI].x
                Δy = nodeList[nodeJ].y - nodeList[nodeI].y
        
                d = radiusList[nodeI] + radiusList[nodeJ]
        
                fx,fy = ForceCalculation(ϵ,d,rc,Δx,Δy)
                cfX[nodeI] += fx
                cfY[nodeI] += fy
            end 
        end
    end

    nothing
end


# ------
# TESTS
# ------

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