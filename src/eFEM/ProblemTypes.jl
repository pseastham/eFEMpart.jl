# Code to be loaded into eFEMpart

############################
###### Boundary Types ######
############################

abstract type AbstractProblem end
abstract type AbstractTimeInfo end
abstract type BoundaryID       end
abstract type BoundaryFunction end

struct DirichletID <: BoundaryID id::Vector{Symbol} end
struct NeumannID <: BoundaryID id::Vector{Symbol} end
struct RobinID <: BoundaryID id::Vector{Symbol} end

struct DirichletFunction <: BoundaryFunction f::Function end
struct NeumannFunction <: BoundaryFunction f::Function end
struct RobinFunction <: BoundaryFunction f::Function end
struct ForcingFunction <: BoundaryFunction f::Function end

#### constructors ####

function Dirichlet(s::Symbol ...)
  return DirichletID([i for i in s])
end

function Neumann(s::Symbol ...)
  return NeumannID([i for i in s])
end

function Robin(s::Symbol ...)
  return RobinID([i for i in s])
end

function Dirichlet(f::Function)
  return DirichletFunction(f)
end

function Neumann(f::Function)
  return NeumannFunction(f)
end

function Robin(f::Function)
  return RobinFunction(f)
end

function Forcing(f::Function)
  return ForcingFunction(f)
end

#### Methods #####

isDirichletID(x) = (typeof(x)==DirichletID ? true : false)
isNeumannID(x)   = (typeof(x)==NeumannID   ? true : false)
isRobinID(x)     = (typeof(x)==RobinID     ? true : false)
isDirichletFunction(x) = (typeof(x)==DirichletFunction ? true : false)
isNeumannFunction(x)   = (typeof(x)==NeumannFunction   ? true : false)
isRobinFunction(x)     = (typeof(x)==RobinFunction     ? true : false)
isForcingFunction(x)   = (typeof(x)==ForcingFunction   ? true : false)

############################
###### Problem Types #######
############################

###############################
######## TimeInfo Types #######
###############################

struct TimeInfo <: AbstractTimeInfo
  TimeDependent::Bool
  tspan::Vector{Float64}
end

##### TimeInfo Constructors #####
# if only given bool, assume it is false, and so tspan is automatically set
TimeInfo(td::Bool) = TimeInfo(td,[0.0,0.0])

mutable struct Problem <: AbstractProblem
  name::String
  tInfo::TimeInfo
  OperatorType::Symbol
  bcid::DictList{Int}
  bcval::DictList{Float64}
end

######  Constructors for Problem  ##########

# Assumes that only 1 dNodes and nNodes and rNodes (maximally ...)
function Problem(mesh::ScalarMesh,Nodes,bcfun,OpType,name,tinfo)
  bcIDNodes::Vector{Vector{Int}} = []
  bcIDNamesSym::Vector{Symbol} = []
  bcValNodes::Vector{Vector{Float64}} = []
  bcValNamesSym::Vector{Symbol} = []

  # derive boundary node ID information
  dNodesSym::Vector{DirichletID} = filter(isDirichletID,Nodes)
  if length(dNodesSym)>0
    dNodes::Vector{Int} = vcat([mesh.bdry[s] for s in dNodesSym[1].id]...)
    push!(bcIDNodes,dNodes)
    push!(bcIDNamesSym,:dNodes)
  end

  nNodesSym::Vector{NeumannID} = filter(isNeumannID,Nodes)
  if length(nNodesSym)>0
    nNodes::Vector{Int} = vcat([mesh.bdry[s] for s in nNodesSym[1].id] ...)
    push!(bcIDNodes,nNodes)
    push!(bcIDNamesSym,:nNodes)
  end

  rNodesSym = filter(isRobinID,Nodes)
  if length(rNodesSym)>0
    rNodes::Vector{Int} = vcat([mesh.bdry[s] for s in rNodesSym[1].id] ...)
    push!(bcIDNodes,rNodes)
    push!(bcIDNamesSym,:rNodes)
  end


  # derive boundary value information
  dNodesFun::Vector{DirichletFunction} = filter(isDirichletFunction,bcfun)
  if length(dNodesFun)>0
    dBCarr::Vector{Float64} = [dNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in dNodes]
    push!(bcValNodes,dBCarr)
    push!(bcValNamesSym,:dBC)
  end

  nNodesFun::Vector{NeumannFunction} = filter(isNeumannFunction,bcfun)
  if length(nNodesFun)>0
    nBCarr::Vector{Float64} = [nNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in nNodes]
    push!(bcValNodes,nBCarr)
    push!(bcValNamesSym,:nBC)
  else
    if :nNodes in bcIDNamesSym
      nBCarr = zeros(Float64,length(nNodes))
      push!(bcValNodes,nBCarr)
      push!(bcValNamesSym,:nBC)
    end
  end

  rNodesFun::Vector{RobinFunction} = filter(isRobinFunction,bcfun)
  if length(rNodesFun)>0
    rBCarr::Vector{Float64} = [rNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in rNodes]
    push!(bcValNodes,rBCarr)
    push!(bcValNamesSym,:rBC)
  end

  fNodesFun::Vector{ForcingFunction} = filter(isForcingFunction,bcfun)
  if length(fNodesFun)>0
    farr::Vector{Float64} = [fNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i=1:length(mesh.xy)]
    push!(bcValNodes,farr)
    push!(bcValNamesSym,:forcing)
  else
    farr = zeros(Float64,length(mesh.xy))
    push!(bcValNodes,farr)
    push!(bcValNamesSym,:forcing)
  end

  # convert all these into appropriate types for problem object
  bcIDNames = Index(bcIDNamesSym)
  bcid  = DictList{Int}(bcIDNodes,bcIDNames)

  bcValNames = Index(bcValNamesSym)
  bcval  = DictList{Float64}(bcValNodes,bcValNames)

  # return FluidProblem object
  return Problem(name,tinfo,OpType,bcid,bcval)
end

function Problem(mesh::ScalarMesh,Nodes,bcfun,name,OpType)
  tinfo = TimeInfo(false)
  return Problem(mesh,Nodes,bcfun,OpType,name,tinfo)
end

# assumes that dUNodes = dVNodes and anything that isn't dirichlet is neumann
function Problem(mesh::FluidMesh,Nodes,bcfun,OperatorType)
  bcIDNodes::Vector{Vector{Int}} = []
  bcIDNamesSym::Vector{Symbol} = []
  bcValNodes::Vector{Vector{Float64}} = []
  bcValNamesSym::Vector{Symbol} = []

  # derive boundary node ID information
  dNodesSym::Vector{DirichletID} = filter(isDirichletID,Nodes)
  if length(dNodesSym)>0
    dUNodes::Vector{Int} = unique(vcat([mesh.bdry[s]
                                  for s in dNodesSym[1].id]...))
    push!(bcIDNodes,dUNodes)
    push!(bcIDNamesSym,:dUNodes)

    if length(dNodesSym)==2
      dVNodes::Vector{Int} = unique(vcat([mesh.bdry[s]
                                    for s in dNodesSym[2].id]...))
    else
      dVNodes = copy(dUNodes)
    end

    push!(bcIDNodes,dVNodes)
    push!(bcIDNamesSym,:dVNodes)
  end

  nNodesSym::Vector{NeumannID} = filter(isNeumannID,Nodes)
  if length(nNodesSym)==1
    nNodes::Vector{Int} = unique(vcat([mesh.bdry[s]
                                 for s in nNodesSym[1].id] ...))
    push!(bcIDNodes,nNodes)
    push!(bcIDNamesSym,:nNodes)
  end

  # derive boundary value information
  dNodesFun::Vector{DirichletFunction} = filter(isDirichletFunction,bcfun)
  if length(dNodesFun)==2
    dUBCarr::Vector{Float64} = [dNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in dUNodes]
    push!(bcValNodes,dUBCarr)
    push!(bcValNamesSym,:dUBC)

    dVBCarr::Vector{Float64} = [dNodesFun[2].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in dVNodes]
    push!(bcValNodes,dVBCarr)
    push!(bcValNamesSym,:dVBC)
  end

  nNodesFun::Vector{NeumannFunction} = filter(isNeumannFunction,bcfun)
  if length(nNodesFun)==1
    nBCarr::Vector{Float64} = [nNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in nNodes]
    push!(bcValNodes,nBCarr)
    push!(bcValNamesSym,:nBC)
  else
    if :nNodes in bcIDNamesSym
      nBCarr = zeros(Float64,length(nNodes))
      push!(bcValNodes,nBCarr)
      push!(bcValNamesSym,:nBC)
    end
  end

  fNodesFun::Vector{ForcingFunction} = filter(isForcingFunction,bcfun)
  if length(fNodesFun)==2
    farrX::Vector{Float64} = [fNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i=1:length(mesh.xy)]
    farrY::Vector{Float64} = [fNodesFun[2].f(node.x,node.y)
                               for node in mesh.xy]
    push!(bcValNodes,farrX)
    push!(bcValNamesSym,:forcingX)
    push!(bcValNodes,farrY)
    push!(bcValNamesSym,:forcingY)
  else
    farrX = zeros(Float64,length(mesh.xy))
    farrY = zeros(Float64,length(mesh.xy))
    push!(bcValNodes,farrX)
    push!(bcValNamesSym,:forcingX)
    push!(bcValNodes,farrY)
    push!(bcValNamesSym,:forcingY)
  end

  # convert all these into appropriate types for problem object
  bcIDNames = Index(bcIDNamesSym)
  bcid  = DictList{Int}(bcIDNodes,bcIDNames)

  bcValNames = Index(bcValNamesSym)
  bcval  = DictList{Float64}(bcValNodes,bcValNames)

  name = "fluid"
  tinfo = TimeInfo(false)

  # return Problem object
  return Problem(name,tinfo,OperatorType,bcid,bcval)
end

####### Problem Methods #######

bc(prob::Problem) = indexlist(prob.bcid)
