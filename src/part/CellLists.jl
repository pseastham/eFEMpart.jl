# code to be loaded into eFEMpart

# This contains functions for generating cell lists for keeping track of
# near-field interaction forces to avoid computational cost of summing
# many unnecessary zeros

# generates list of cells,
# where rm is the width of each cell block
function CellGenerate(particleList::Vector{Particle},
                      square::Vector{Float64},rm::Float64)
  X0=square[1]; X1=square[2]; Y0=square[3]; Y1=square[4]
  Nx = Int(ceil((X1 - X0)/rm))
  Ny = Int(ceil((Y1 - Y0)/rm))
  Nparticles = length(particleList)

  CellList = Array{Cell}(undef,Nx*Ny)

  k=1
  for i=0:(Nx-1), j=1:Ny
    # compute square
    x0 = X0 + i*rm; x1 = X0 + (i+1)*rm; y0 = Y1-j*rm; y1 = Y1-(j-1)*rm

    # compute neighbors
    neighbors = []
    # left
    if i>0 push!(neighbors,(i-1)*Nx+j) end
    # top
    if j>1 push!(neighbors,i*Nx+j-1) end
    # right
    if i<(Nx-1) push!(neighbors,(i+1)*Nx+j) end
    # bottom
    if j<Ny push!(neighbors,i*Nx+j+1) end
    # top-right
    if (j>1 && i<(Nx-1)) push!(neighbors,(i+1)*Nx+j-1) end
    # top-left
    if (j>1 && i>0) push!(neighbors,(i-1)*Nx+j-1) end
    # bottom-right
    if (j<Ny && i<(Nx-1)) push!(neighbors,(i+1)*Nx+j+1) end
    # bottom-left
    if (j<Ny && i>0) push!(neighbors,(i-1)*Nx+j+1) end

    # generate cell
    CellList[k] = Cell([],[x0,x1,y0,y1],neighbors)
    k+=1
  end

  # assign particles
  for i=1:Nparticles
    point = [particleList[i].xpos,particleList[i].ypos]
    for l=1:Nx*Ny
      if isInsideRect(CellList[l].square,point)
        push!(CellList[l].ParticleList,i)
        break
      end
    end
  end

  return CellList
end

# generates cell list for a mesh,
# useful when interpolating the FEM solution
# again rm is the length of the cell block
function CellGenerate(mesh::AbstractMesh,rm::Float64)
  Nnodes = length(mesh.xy)

  # generate particle List for FEM nodes
  pList = [Particle(mesh.xy[i].x,mesh.xy[i].y) for i=1:Nnodes]

  # generate square
  x0 = 0.0; x1=0.0; y0=0.0; y1=0.0

  # generate cell list
  for i=1:Nnodes
    x = mesh.xy[i].x; y = mesh.xy[i].y
    if x < x0
      x0 = x
    elseif x > x1
      x1 = x
    end
    if y < y0
      y0 = y
    elseif y > y1
      y1 = y
    end
  end

  xlen = x1-x0
  ylen = y1-y0

  square = [x0-0.1*xlen,x1+0.1*xlen,y0-0.1*ylen,y1+0.1*ylen]
  # return cell list
  femCL = CellGenerate(pList,square,rm)

  return femCL
end

# determines whether or not point is inside a rectangle
function isInsideRect(rect::Vector{Float64},point::Vector{Float64})
  x=point[1]; y=point[2]
  x0=rect[1]; x1=rect[2]; y0=rect[3]; y1=rect[4]

  if x<x1 && x>x0 && y<y1 && y>y0
    return true
  else
    return false
  end
end
