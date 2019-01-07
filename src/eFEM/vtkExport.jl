  # File to be put into eFEM module

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
  New scalar + vector solutions -- not used. outdated.
"""
function vtksave2(mesh::T,sd::ScalarData,sn::ScalarNames,vd::VectorData,vn::VectorNames,path::Path;time=0.0) where T<:AbstractMesh
  Nscalars = length(sd.scalarArr)
  Nvectors = length(vd.vectorArr)

  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  # convert pressure from Q1 to Q2
  pressureID = 0
  for ti=1:Nscalars
    if sn.nameArr[ti] == "pressure"
      pressureID = ti
      break
    end
  end
  if pressureID != 0 && typeof(mesh)==FluidMesh
      # compute pressure at Q2 mesh
      xyArr = NodesToArray(mesh.xy)
      xypArr = NodesToArray(mesh.xyp)
      cmpArr = ElementsToArray(mesh.cmp)
      newp = meshTransform(xypArr,cmpArr,sd.scalarArr[pressureID],xyArr)
      sd.scalarArr[pressureID] = newp
  end

  newscalar = zeros(Float64,length(mesh.xy),Nscalars)
  for i=1:length(mesh.xy), j=1:Nscalars
      newscalar[i,j] = sd.scalarArr[j][i]
  end

  newvector = zeros(Float64,length(mesh.xy),2,Nvectors)
  for i=1:length(mesh.xy), j=1:2, k=1:Nvectors
    newvector[i,j,k] = vd.vectorArr[k][j][i]
  end

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`);
  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,newscalar,sn.nameArr,newvector,vn.nameArr,file_name;time=time)
end
"""
  New scalar + vector solutions -- ONLY Q1 MESH -- use this one!!
"""
function vtksave(mesh::T,sd::ScalarData,sn::ScalarNames,vd::VectorData,vn::VectorNames,path::Path;time=0.0) where T<:AbstractMesh
  if true #mesh is Q2!!!! (or fluid)
  end
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

  if typeof(mesh) == FluidMesh
    newscalar = zeros(Float64,size(xyp,1),Nscalars)
    for j=1:Nscalars
      if sn.nameArr[j] != "pressure"
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
  elseif mesh.order == :Linear
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
  else
    error("this wasn't set up!") 
  end

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`);
  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xyp,cmp,newscalar,sn.nameArr,newvector,vn.nameArr,file_name;time=time)
end

# ===================================================================== #
# ===================================================================== #
# ===================================================================== #
# ===================================================================== #

"""
  VTK output for standard array solution
"""
function vtksave(mesh::T,sol::Vector{S},name::String,path::Path;time=0.0) where
                 {T<:AbstractMesh,S<:Real}
  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`);
  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,sol,name,file_name;time=time)
end
"""
  VTK export for multiple array solutions
"""
function vtksave(mesh::T,sol::Vector{Vector{S}},names::Vector{String},path::Path;time=0.0) where 
                {T<:AbstractMesh,S<:Real}
  Nsols = length(sol)


  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  newsol = zeros(length(mesh.xy),Nsols)
  for i=1:length(mesh.xy)
    for j=1:Nsols
      newsol[i,j] = sol[j][i]
    end
  end

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`);
  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,newsol,names,file_name;time=time)

end
"""
  VTK output for my ScalarSolution type
"""
function vtksave(mesh::ScalarMesh,prob::Problem,param::T,sol::ScalarSolution,
                 path::Path) where T<:AbstractParameter
  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`)

  # potential velocity field
  if typeof(param) <: AbstractAdvDiff
    vName = "velocity"; vValtemp  = hcat(param.u,param.v)
    vVal = reshape(vValtemp,size(vValtemp)[1],size(vValtemp)[2],1)
  else
    vName = []; vVal  = []
  end

  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,sol.u,prob.name,vVal,vName,file_name)
end
"""
  VTK output for my FluidSolution type
"""
function vtksave(mesh::FluidMesh,prob,param::T,sol::FluidSolution,
                 path::Path) where T<:AbstractParameter
  run(`mkdir -p $(path.folder)`)

  file_name = string(path.folder,"/",path.file_name,".vtk")
  fluidVTK(mesh,sol;fileName=file_name,sName=["other"],sVal=[])
end
"""
  VTK output without path
"""
function vtksave(mesh::M,prob::AbstractProblem,param::T,sol::S) where 
                 {M<:AbstractMesh,T<:AbstractParameter,S<:AbstractSolution}
  path = Path("solutions/",string("FEM",prob.name))
  return vtksave(mesh,prob,param,sol,path)
end
"""
  VTK output without parameter (insert dummy parameter for :Poisson2D)
"""
function vtksave(mesh::M,prob::P,sol::S) where
                 {M<:AbstractMesh,P<:AbstractProblem,S<:AbstractSolution}
  param = PoissonParam(1.0)
  return vtksave(mesh,prob,param,sol)
end
"""
  VTK output for multiple solutions all using same parameter set

  In the future could possibly reduce code by tying this and the next 
  one together
"""
function vtksave(mesh::M,probArr::Vector{Problem},param::T,
                 solArr::Vector{ScalarSolution},path::Path) where 
                 {M<:AbstractMesh,T<:AbstractAdvDiff}        
  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`)

  vName = "velocity"; vVal = hcat(param.u,param.v)

  if length(probArr>3); error("vtksave not built to handle other than 3 inputs.") end
  fNames = hcat(probArr[1].name,probArr[2].name,probArr[3].name)
  fvals  = hcat(solArr[1].u,solArr[2].u,solArr[3].u)

  file_name = string(path.folder,"/",path.file_name,".vtk")
  return TWOD_UNSTRUCTURED_TO_VTK(xy,cm,fvals,fNames,vVal,vName,file_name)
end
"""
  VTK output for multiple solutions all using same parameter set

  In the future could possibly reduce code by tying this and the previous 
  one together
"""
function vtksave(mesh::M,probArr::Vector{Problem},param::T,
                 solArr::Vector{Vector{R}},path::Path) where 
                 {M<:AbstractMesh,T<:AbstractAdvDiff,R<:Real}   
  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`)

  vName = "velocity"; vValtemp  = hcat(param.u,param.v)
  vVal = reshape(vValtemp,size(vValtemp)[1],size(vValtemp)[2],1)

  if length(probArr>3); error("vtksave not built to handle other than 5 inputs.") end
  fNames = hcat(probArr[1].name,probArr[2].name,probArr[3].name,
                probArr[4].name,probArr[5].name)
  fvals  = hcat(solArr[1],solArr[2],solArr[3],solArr[4],solArr[5])

  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,fvals,fNames,vVal,vName,file_name)
end
"""
  VTK output for Darcy's equation
"""
function vtksave(mesh::ScalarMesh,prob::Problem,
                 param::T,sol::FluidSolution,
                 path::Path) where T<:AbstractDarcy
  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  # creates folder if it doesn't already exist
  run(`mkdir -p $(path.folder)`)

  fNames = prob.name; fvals = sol.p
  vName = "velocity"; vValtemp  = hcat(sol.u,sol.v)
  vVal = reshape(vValtemp,size(vValtemp)[1],size(vValtemp)[2],1)

  file_name = string(path.folder,"/",path.file_name,".vtk")
  TWOD_UNSTRUCTURED_TO_VTK(xy,cm,fvals,fNames,vVal,vName,file_name)
end
"""
  VTK output for plain mesh
"""
function vtksave(mesh::T,path::Path) where T<:AbstractMesh
    # get xy & cm into correct form
    xy = NodesToArray(mesh.xy)
    cm = ElementsToArray(mesh.cm)
  
    # creates folder if it doesn't already exist
    run(`mkdir -p $(path.folder)`)
  
    fNames = ""; fvals = zeros(Float64,length(mesh.xy))
    vName = ""; vVal = []
  
    file_name = string(path.folder,"/",path.file_name,".vtk")
    TWOD_UNSTRUCTURED_TO_VTK(xy,cm,fvals,fNames,vVal,vName,file_name)
end
"""
  VTK output for plain mesh -- name as String, not Path object
"""
function vtksave(mesh::T,pathString::String) where T<:AbstractMesh
  path = Path(pathString)
  return vtksave(mesh,path)
end
"""
  VTK output for plain mesh -- default name
"""
function vtksave(mesh::T) where T<:AbstractMesh
    return vtksave(mesh,"mesh")
end
"""
  VTK output for particles
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
  fluidVTK(u,v,p,xy,cm,xyp,cmp;
           folder='solutions',fileName='fluidSolution',
           sName='other',sVal=[])

Outputs solution (u,v,p) to VTK file.

Currently can accomidate additional scalar field using sName and sVal optional fields

LATER: include option to compute vorticity and stress fields
"""
function fluidVTK(mesh::FluidMesh,sol=FluidSolution;
                  fileName="solutions/fluidSolution.vtk",
                  sName=["other"],sVal=[])
  # get xy & cm into correct form
  xy = NodesToArray(mesh.xy);    xyp = NodesToArray(mesh.xyp)
  cm = ElementsToArray(mesh.cm); cmp = ElementsToArray(mesh.cmp)

  u = sol.u; v = sol.v; p = sol.p

  q1u,q1v = lowerUV(xy,cm,u,v)
  vValtemp = hcat(q1u,q1v)
  vel = reshape(vValtemp,size(vValtemp)[1],size(vValtemp)[2],1)

  fName = ["pressure"]
  fVal  = p

  if typeof(sName)==String
    q1s     = lowerScalar(xy,cm,sVal)

    fName = hcat(fName,sName)
    fVal  = hcat(fVal,q1s)
  elseif typeof(sName) == Array{String,1}
    if size(sVal,1) != 0
      if length(sName) == 1
        q1s   = lowerScalar(xy,cm,sVal[i])

        fName = hcat(fName,sName[i])
        fVal  = hcat(fVal,q1s)
      else
        for i=1:length(sName)
          q1s   = lowerScalar(xy,cm,sVal[:,i])

          fName = hcat(fName,sName[i])
          fVal  = hcat(fVal,q1s)
        end
      end
    end
  end

  vecname = "velocity"
  TWOD_UNSTRUCTURED_TO_VTK(xyp,cmp,fVal,fName,vel,vecname,fileName)
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
