# computes cohesion forces

include("CellLists.jl")  # loads in Point2D type

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

function computeCohesion!(cfX::Vector{T},cfY::Vector{T},nodeList::Vector{Point2D{T}},radiusList::Vector{T},rc::T,ϵ::T) where T<:Real
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
function computeCohesion_CL!(cfX::Vector{T},cfY::Vector{T},nodeList::Vector{Point2D{T}},radiusList::Vector{T},
                                rc::T,ϵ::T,cl::CellList) where T<:Real
    # zero-out cohesion
    fill!(cfX,zero(T))
    fill!(cfY,zero(T))

    # go over cell lists
    for cellInd = 1:length(cl.cells)
        # compute cell-cell (same cell) forces
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