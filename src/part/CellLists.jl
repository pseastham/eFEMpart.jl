# Generates and updates cell list object, and generates mapping for 
# finite element mesh elements to a corresponding cell list cells

include("ParticleTypes.jl")  # to load in Point2D type

# -----------------
# TYPE DEFINITIONS 
# -----------------

mutable struct Cell
  nodeList::Vector{Int}           # list of indices of particles that are within cell
  square::Vector{Float64}         # bounding square that defined geometry and position of cell
  neighborList::Vector{Int}       # list of cell neighbors, used for fast access
end

mutable struct CellList
  cells::Vector{Cell}             # array of cells -- this is main component
  bounds::Vector{Float64}         # total bounds of cell list in [x0, x1, y0, y1] format
  sideLength::Float64             # length of one cell inside cell list
end

# ------------
# CONSTRUCTOR 
# ------------

"""
  generateCellList(nodeList,square,L)

  Generates list of cells, where L is the width of each cell block

INPUT
nodeList
TotalBounds
L

OUTPUT
cl
"""
function generateCellList(nodeList::Vector{Point2D},TotalBounds::Vector{Float64},L::Float64)
  X0=TotalBounds[1]; X1=TotalBounds[2]; Y0=TotalBounds[3]; Y1=TotalBounds[4]
  Nx = Int(ceil((X1 - X0)/L))
  Ny = Int(ceil((Y1 - Y0)/L))
  nNodes = length(nodeList)

  neighborMat = zeros(Int,8)
  cl = CellList(Array{Cell}(undef,Nx*Ny),TotalBounds,L)

  tempZeroArr = Array{Int}(undef,0)

  # compute neighbors of cells
  k=1
  for i=0:(Nx-1), j=1:Ny
    matInd = 1

    # zero-out neighborMat
    fill!(neighborMat,zero(Int))

    # compute bounding square
    x0 = X0 + i*L; x1 = X0 + (i+1)*L; y0 = Y1-j*L; y1 = Y1-(j-1)*L

    # compute neighbors
    # left
    if i>0 neighborMat[matInd] = (i-1)*Nx+j; matInd+=1 end
    # top
    if j>1 neighborMat[matInd] = i*Nx+j-1; matInd+=1 end
    # right
    if i<(Nx-1) neighborMat[matInd] = (i+1)*Nx+j; matInd+=1 end
    # bottom
    if j<Ny neighborMat[matInd] = i*Nx+j+1; matInd+=1 end
    # top-right
    if (j>1 && i<(Nx-1)) neighborMat[matInd] = (i+1)*Nx+j-1; matInd+=1 end
    # top-left
    if (j>1 && i>0) neighborMat[matInd] =(i-1)*Nx+j-1; matInd+=1 end
    # bottom-right
    if (j<Ny && i<(Nx-1)) neighborMat[matInd] = (i+1)*Nx+j+1; matInd+=1 end
    # bottom-left
    if (j<Ny && i>0) neighborMat[matInd] = (i-1)*Nx+j+1; matInd+=1 end

    # generate cell
    cl.cells[k] = Cell(copy(tempZeroArr),[x0,x1,y0,y1],neighborMat[1:(matInd-1)])
    k+=1
  end

  # assign nodes to respective cells
  for nodeInd=1:nNodes
    for cellInd=1:Nx*Ny
      if isInsideRect(cl.cells[cellInd].square,nodeList[nodeInd])
        push!(cl.cells[cellInd].nodeList,nodeInd)
        break
      end
    end
  end

  return cl
end

"""

"""
function updateCellList!(cl::CellList,nodeList::Vector{Point2D})
  X0=cl.bounds[1]; X1=cl.bounds[2]; Y0=cl.bounds[3]; Y1=cl.bounds[4]
  Nx = Int(ceil((X1 - X0)/cl.sideLength))
  Ny = Int(ceil((Y1 - Y0)/cl.sideLength))
  nNodes = length(nodeList)

  # 'zero-out' cell node lists
  for cellInd=1:Nx*Ny
    cl.cells[cellInd].nodeList = Array{Int}(undef,0)
  end

  # assign nodes to respective cells
  for nodeInd=1:nNodes
    for cellInd=1:Nx*Ny
      if isInsideRect(cl.cells[cellInd].square,nodeList[nodeInd])
        push!(cl.cells[cellInd].nodeList,nodeInd)
        break
      end
    end
  end

  nothing
end

"""

updates cell list with assumption that most nodes are in same cell as before, or have gone to neighbors

should be much faster, but isn't......

DO NOT USE
"""
function updateCellList_SMART!(cl::CellList,nodeList::Vector{Point2D})
  X0=cl.bounds[1]; X1=cl.bounds[2]; Y0=cl.bounds[3]; Y1=cl.bounds[4]
  Nx = Int(ceil((X1 - X0)/cl.sideLength))
  Ny = Int(ceil((Y1 - Y0)/cl.sideLength))
  nNodes = length(nodeList)

  nodesToFind = Array{Int}(undef,0)

  # assign nodes to respective cells
  for cellInd=1:Nx*Ny
    numDeleted = 0
    for ti=1:length(cl.cells[cellInd].nodeList)
      found = false
      if !(isInsideRect(cl.cells[cellInd].square,nodeList[cl.cells[cellInd].nodeList[ti-numDeleted]]))
        # search neighbors of cell for node
        for nbrs in cl.cells[cellInd].neighborList
          if isInsideRect(cl.cells[nbrs].square,nodeList[cl.cells[cellInd].nodeList[ti-numDeleted]])
            push!(cl.cells[nbrs].nodeList,cl.cells[cellInd].nodeList[ti-numDeleted])
            found = true
            break
          end
        end

        if !(found)
          push!(nodesToFind,cl.cells[cellInd].nodeList[ti-numDeleted])
        end
        # delete node from cell nodelist
        #println("node deleted")
        deleteat!(cl.cells[cellInd].nodeList,ti-numDeleted)
        numDeleted += 1
      end
    end
  end

  # for nodes not found, search through all elements
  #println(length(nodesToFind))
  for nodeInd in nodesToFind
    for cellInd=1:Nx*Ny
      if isInsideRect(cl.cells[cellInd].square,nodeList[nodeInd])
        push!(cl.cells[cellInd].nodeList,nodeInd)
        break
      end
    end
  end

  nothing
end

# generates array of arrays for a mesh, useful when interpolating the FEM solution
# L is the length of a single cell
#
# index of femCLmap corresponds to cell index
# values of femCLmap[jj] correspond to index of elements inside cell jj
function femGenerateMap(mesh::M,totalBounds::Vector{Float64},L::Float64) where M#<:AbstractMesh
  Nnodes = length(mesh.xy)

  # generate particle List for FEM nodes
  nodeList = [Point2D(mesh.xy[i].x,mesh.xy[i].y) for i=1:Nnodes]

  tempCL = generateCellList(nodeList,totalBounds,L)
  femCLmap = Array{Array{Int}}(undef,length(tempCL.cells))

  # initialize element array
  elementArray = Array{Int}(undef,4)

  # convert all point indices to corresponding element indices
  for cellInd = 1:length(tempCL.cells)
    tempArr = Array{Int}(undef,0)

    #for each node in that cell
    for nodeInd in tempCL.cells[cellInd].nodeList
      fill!(elementArray,zero(Int))

      # map node index to element
      getElementFromNode!(elementArray,mesh,nodeInd)

      for ti=1:4
        if elementArray[ti] != 0
          push!(tempArr,elementArray[ti])
        end
      end
    end

    # assign elements to femCL as 'nodes'
    femCLmap[cellInd] = unique(tempArr)
  end

  return femCLmap
end

# determines whether or not point is inside a rectangle
# also includes being on bottom or left boundary
function isInsideRect(rect::Vector{Float64},point::Point2D)
  x0=rect[1]; x1=rect[2]; y0=rect[3]; y1=rect[4]

  if point.x<x1 && point.x>=x0 && point.y<y1 && point.y>=y0
    return true
  else
    return false
  end
end

function getElementFromNode!(elementArray::Vector{Int},mesh::M,nodeInd::Int) where M#<:AbstractMesh
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