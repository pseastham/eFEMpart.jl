# generates array of arrays for a mesh, useful when interpolating the FEM solution
# L is the length of a single cell
#
# index of femCLmap corresponds to cell index
# values of femCLmap[jj] correspond to index of elements inside cell jj
function generate_femap(mesh::M,totalBounds::Vector{Float64},L::Float64) where M<:AbstractMesh
    Nnodes = length(mesh.xy)
    numNodesPerElm = (mesh.order == :Linear ? 4 : 9)

    # generate particle List for FEM nodes
    particleList = [Particle2D(Point2D(mesh.xy[i].x,mesh.xy[i].y),0.0,0) for i=1:Nnodes]

    tempCL = StokesParticles.generate_cell_list(particleList,totalBounds,L)
    femap = Array{Array{Int}}(undef,length(tempCL.cells))

    # initialize element array
    elementArray = Array{Int}(undef,numNodesPerElm)

    # convert all point indices to corresponding element indices
    for cellInd = 1:length(tempCL.cells)
        tempArr = Array{Int}(undef,0)

        #for each node in that cell
        for nodeInd in tempCL.cells[cellInd].particleIDList
            fill!(elementArray,zero(Int))

            # map node index to element
            getElementFromNode!(elementArray,mesh,nodeInd)

            for ti=1:numNodesPerElm
                if elementArray[ti] != 0
                    push!(tempArr,elementArray[ti])
                end
            end
        end

        # assign elements to femCL as 'nodes'
        femap[cellInd] = unique(tempArr)
    end

    return femap
end 

"""
    interpFGT!(v,pList,mesh,h,ε,x,y,q)

Takes in particle position and FEM mesh to interpolate values.
Also updates particle possitoin arrays in-place to reduce memory consumption

POSSIBLE IMPROVEMENTS:
1. add another function that updates the target array only once the first time
       this code is run. Maybe some sort of interp object that does this for me

INPUTS:
v:     [Vector{Float64}, M] interpolated values at targets
pList: [ParticleList] Array of particles, acts as sources
mesh:  [Mesh] Finite Element mesh, acts as targets
h:     [Float64] bandwidth, h=sqrt(2)*stddev if standard deviation is known
ε:     [Float64] maximum absolute error
x:     [Vector{Float64}, d*N] list of sources -- updated in-place within function
y:     [Vector{Float64}, d*M] list of targets -- updated in-place within function
q:     [Vector{Float64}, N] weights on sources -- updated in-place within function
W:     [Int] number of weight functions (e.g. if W=2, then length(q) = 2N)
"""
function interpFGT!(v::Array{Float64},pList,mesh,h::Float64,ε::Float64,
                    x::Array{Float64},y::Array{Float64},q::Array{Float64})
    N = length(pList)       # sources -- needed for fgt function
    M = length(mesh.xy)     # targets -- needed for fgt function
    d = 2                   # dimension of data -- needed for fgt function
    W = 1                   # number of different weights -- needed for fgt function

    # update source array
    for ti=1:N
        x[2*ti-1] = pList[ti].x
        x[2*ti]   = pList[ti].y
    end

    # update target array -- really should only have to do this once! room for improvement
    for ti=1:M
        y[2*ti-1] = mesh.xy[ti].x       # xval of source
        y[2*ti]   = mesh.xy[ti].y
    end

    # run interpolation routine
    fgt!(v,d,M,N,h,ε,x,y,q,W)

    # return nothing -- updates in-place
    nothing
end

function getElementFromNode!(elementArray::Vector{Int},mesh::M,nodeInd::Int) where M<:AbstractMesh
    nElm = length(mesh.cm)
  
    tk = 1
    for elmInd=1:nElm
      if nodeInd in mesh.cm[elmInd].NodeList
        elementArray[tk] = elmInd
        tk += 1
      end
    end
  
    nothing
  end