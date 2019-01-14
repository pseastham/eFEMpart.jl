  # File to be put into eFEMpart module

struct Path
  folder::String
  file_name::String
end

struct ScalarData 
  scalarArr::Vector{Vector{Float64}}
end

struct ScalarNames
  nameArr::Vector{String}
end

struct VectorData
  vectorArr::Vector{Vector{Vector{Float64}}}
end

struct VectorNames
  nameArr::Vector{String}
end

# -----------
# constructors
# -----------

# single entry
ScalarData(u::Vector{Float64}) = ScalarData([u])
ScalarNames(u::String) = ScalarNames([u])
VectorData(u::Vector{Vector{Float64}}) = VectorData([u])
VectorNames(u::String) = VectorNames([u])

# multiple entries
ScalarData(u::Vector{Float64}...) = ScalarData([u[i] for i=1:length(u)])
ScalarNames(u::String...) = ScalarNames([u[i] for i=1:length(u)])
VectorData(u::Vector{Vector{Float64}}...) = VectorData([u[i] for i=1:length(u)])
VectorNames(u::String...) = VectorNames([u[i] for i=1:length(u)])

Path(file_name::String) = Path("solutions",file_name)

##########################
#### vtsave functions ####
##########################

"""
  New scalar solutions
"""
function vtksave(mesh::T,sd::ScalarData,sn::ScalarNames,path::Path;time=0.0) where T<:AbstractMesh
  Nscalars = length(sd.scalarArr)

  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  newsol = zeros(Float64,length(mesh.xy),Nscalars)
  for i=1:length(mesh.xy), j=1:Nscalars
    newsol[i,j] = sd.scalarArr[j][i]
  end

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`);
  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,newsol,sn.nameArr,file_name;time=time)
end
"""
  New vector solutions
"""
function vtksave(mesh::T,vd::VectorData,vn::VectorNames,path::Path;time=0.0) where T<:AbstractMesh
  Nvectors = length(vd.vectorArr)

  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  newsol = zeros(Float64,length(mesh.xy),2,Nvectors)
  for i=1:length(mesh.xy), j=1:2, k=1:Nvectors
    newsol[i,j,k] = vd.vectorArr[k][j][i]
  end

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`);
  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,[],[],newsol,vn.nameArr,file_name;time=time)
end
"""
  New scalar + vector solutions -- exports linear element visualization
"""
function vtksave(mesh::T,sd::ScalarData,sn::ScalarNames,vd::VectorData,vn::VectorNames,path::Path;time=0.0) where T<:AbstractMesh
  #if true #mesh is Q2!!!! (or fluid)
  #end
  Nscalars = length(sd.scalarArr)
  Nvectors = length(vd.vectorArr)

  # get xy & cm into correct form
  if typeof(mesh) == FluidMesh
    xy = NodesToArray(mesh.xy)
    cm = ElementsToArray(mesh.cm)
    xyp = NodesToArray(mesh.xyp)
    cmp = ElementsToArray(mesh.cmp)
  else
    xy = NodesToArray(mesh.xy)
    cm = ElementsToArray(mesh.cm)
    xyp = NodesToArray(mesh.xy)
    cmp = ElementsToArray(mesh.cm)
  end

  # if fluid mesh, account for pressure
  if typeof(mesh) == FluidMesh
    newscalar = zeros(Float64,size(xyp,1),Nscalars)
    for j=1:Nscalars
      if (sn.nameArr[j] != "pressure") && (sn.nameArr[j] != "head")
        newscalar[:,j] = lowerScalar(xy,cm,sd.scalarArr[j])
      else
        newscalar[:,j] = sd.scalarArr[j]
      end
    end

    newvector = zeros(Float64,size(xyp,1),2,Nvectors)
    for k=1:Nvectors
      uL,vL = lowerUV(xy,cm,vd.vectorArr[k][1],vd.vectorArr[k][2])
      for ti=1:size(xyp,1)
        newvector[ti,1,k] = uL[ti]
        newvector[ti,2,k] = vL[ti]
      end
    end
  else
    newscalar = zeros(Float64,size(xyp,1),Nscalars)
    for j=1:Nscalars
      newscalar[:,j] = sd.scalarArr[j]
    end

    newvector = zeros(Float64,size(xyp,1),2,Nvectors)
    for k=1:Nvectors
      for ti=1:size(xyp,1)
        newvector[ti,1,k] = vd.vectorArr[k][1][ti]
        newvector[ti,2,k] = vd.vectorArr[k][2][ti]
      end
    end  
  end

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`);
  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xyp,cmp,newscalar,sn.nameArr,newvector,vn.nameArr,file_name;time=time)
end

"""
  vtksave for particle output
"""
function vtksave(pList::Vector,radius::Float64,fn::String,timeval::Float64)
  Nparticles = length(pList)
  xy = zeros(Float64,Nparticles,2)
  for i=1:Nparticles
    xy[i,1] = pList[i].xpos
    xy[i,2] = pList[i].ypos
  end

  PARTICLES_TO_VTK(xy,radius,fn;time=timeval)
end

"""
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,fval,fname,vecVal,vecName,file_name[;time=0.0])

  Outputs visualization file in vtk legacy format from FEM data 

Arguments:
`xy::Matrix{Float64}`
`cm::Matrix{Float64}`
`fval`
`fname`
`vecVal`
`vecName`
`file_name::String`
`time::Float64`
"""
function TWOD_UNSTRUCTURED_TO_VTK(xy::Matrix{Float64},cm::Matrix{Int},
                                  fval,fname,vecVal,vecName,
                                  file_name::String;time=0.0)
  # prepare
  fvalPresent   = false
  vecValPresent = false
  Npts = size(xy,1)
  Nels = size(cm,1)
  NPE = size(cm,2)      # nodes per element
  if NPE == 4
    NPEcode = 9
  elseif NPE == 9
    NPE = 8
    NPEcode = 23
  else
    error("Invalid cell type. Must be linear or biquadratic quadrilateral")
  end

  # check that fnames is scalar, and turns it into array of scalars
  if typeof(fname) == String
    fname = [fname]
  end

  # check that fvals or vecVals are nonempty
  if size(fval,1) == Npts
    Numfvals    = size(fval,2)
    fvalPresent = true
  end
  if size(vecVal,1) == Npts
    NumVecVals = size(vecVal,3)

    vecValPresent = true
  end

  # open file to write
  f = open(file_name,"w")

  # header
  println(f,"# vtk DataFile Version 3.0")
  println(f,"")
  println(f,"ASCII")

  printBlankLine(f)

  # print domain points
  println(f,"DATASET UNSTRUCTURED_GRID")

  # print time
  println(f,"FIELD FieldData 1")
  println(f,"TIME 1 1 double")
  println(f,time)

  # print points
  println(f,"POINTS ",Npts," float")
  for i=1:Npts
    println(f,xy[i,1]," ",xy[i,2]," ",0.0)
  end

  printBlankLine(f)

  # print cell information
  println(f,"CELLS ",Nels," ",Nels*(NPE+1))
  for i=1:Nels
    print(f,NPE)
    for elNode = 1:NPE
      print(f," ",cm[i,elNode]-1)
    end
    printBlankLine(f)
  end

  printBlankLine(f)

  # print cell type information
  println(f,"CELL_TYPES ",Nels)
  for i=1:Nels
    println(f,NPEcode)
  end

  printBlankLine(f)

  println(f,"POINT_DATA ",Npts)
  # print scalar values
  if fvalPresent
    for j=1:Numfvals
      println(f,"SCALARS ",fname[j]," double")
      println(f,"LOOKUP_TABLE default")
      for i=1:Npts
        println(f,fval[i,j])
      end
    end
  end

  # print vector values -- 2D vectors
  if vecValPresent
    for j=1:NumVecVals
      println(f,"VECTORS ",vecName[j]," double")
      for i=1:Npts
        println(f,vecVal[i,1,j]," ",vecVal[i,2,j]," ",0.0)
      end
    end
  end

  printBlankLine(f)

	# close file
	close(f)
end
function TWOD_UNSTRUCTURED_TO_VTK(xy::Matrix{Float64},cm::Matrix{Int},
                                  fval,fname,file_name::String;time=0.0)
  vName = []; vVal  = []
  return TWOD_UNSTRUCTURED_TO_VTK(xy,cm,fval,fname,vVal,vName,file_name;time=time)
end

"""
    PARTICLES_TO_VTK(xy,radius,file_name)

Saves particle data in legacy *.vtk file format

# Arguments
- `xy::Matrix{Float64}`: Nx2 matrix of particle positions
- `radius::Float64`: radius of particles
- `file_name::String`: file name for export
"""
function PARTICLES_TO_VTK(xy::Matrix{Float64},radius::Float64,file_name::String;time=0.0)
  # prepare
  Npts = size(xy,1)

  # open file to write
  fnvtk = string(file_name,".vtk")
  f = open(fnvtk,"w")

  # header
  println(f,"# vtk DataFile Version 3.0")
  printBlankLine(f)
  println(f,"ASCII")
  printBlankLine(f)

  # print domain points
  println(f,"DATASET UNSTRUCTURED_GRID")

  # print time
  println(f,"FIELD FieldData 1")
  println(f,"TIME 1 1 double")
  println(f,time)

  # print points  
  println(f,"POINTS ",Npts," float")
  for i=1:Npts
    println(f,xy[i,1]," ",xy[i,2]," ",0.0)
  end

  printBlankLine(f)

  # print cell type information
  println(f,"CELL_TYPES ",Npts)
  for i=1:Npts
    println(f,"1")
  end

  printBlankLine(f)

  println(f,"POINT_DATA ",Npts)
  # print scalar values
  println(f,"SCALARS ","radius"," double")
  println(f,"LOOKUP_TABLE default")
  for i=1:Npts
    println(f,radius)
  end

	# close file
	close(f)
end
