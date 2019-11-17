# Code to be loaded into eFEMpart

############################
###### Boundary Types ######
############################

abstract type AbstractProblem  end
abstract type BoundaryID       end
abstract type BoundaryFunction end
abstract type BoundaryArray end

struct DirichletID <: BoundaryID id::Vector{Symbol} end
struct NeumannID   <: BoundaryID id::Vector{Symbol} end
struct RobinID     <: BoundaryID id::Vector{Symbol} end

struct DirichletFunction <: BoundaryFunction f::Function end
struct NeumannFunction   <: BoundaryFunction f::Function end
struct RobinFunction     <: BoundaryFunction f::Function end
struct ForcingFunction   <: BoundaryFunction f::Function end

struct DirichletArray <: BoundaryArray v::Vector end
struct NeumannArray   <: BoundaryArray v::Vector end
struct RobinArray     <: BoundaryArray v::Vector end
struct ForcingArray   <: BoundaryArray v::Vector end

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
function Dirichlet(v::Vector)
  return DirichletArray(v)
end

function Neumann(f::Function)
  return NeumannFunction(f)
end
function Neumann(v::Vector)
  return NeumannArray(v)
end

function Robin(f::Function)
  return RobinFunction(f)
end
function Robin(v::Vector)
  return RobinArray(v)
end

function Forcing(f::Function)
  return ForcingFunction(f)
end
function Forcing(v::Vector)
  return ForcingArray(v)
end

#### Methods #####

isDirichletID(x)       = (typeof(x)==DirichletID ? true : false)
isNeumannID(x)         = (typeof(x)==NeumannID   ? true : false)
isRobinID(x)           = (typeof(x)==RobinID     ? true : false)

isDirichletFunction(x) = (typeof(x)==DirichletFunction ? true : false)
isNeumannFunction(x)   = (typeof(x)==NeumannFunction   ? true : false)
isRobinFunction(x)     = (typeof(x)==RobinFunction     ? true : false)
isForcingFunction(x)   = (typeof(x)==ForcingFunction   ? true : false)

isDirichletArray(x) = (typeof(x)==DirichletArray ? true : false)
isNeumannArray(x)   = (typeof(x)==NeumannArray   ? true : false)
isRobinArray(x)     = (typeof(x)==RobinArray     ? true : false)
isForcingArray(x)   = (typeof(x)==ForcingArray   ? true : false)

############################
###### Problem Types #######
############################

mutable struct Problem <: AbstractProblem
  OperatorType::Symbol
  bcid::DictList{Int}
  bcval::DictList{Float64}
end

######  Constructors for Problem  ##########

# Assumes that only 1 dNodes and nNodes and rNodes (maximally ...)
function Problem(mesh::ScalarMesh,Nodes,bcfun,OpType)
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
  # dirichlet
  dNodesFun::Vector{DirichletFunction} = filter(isDirichletFunction,bcfun)
  if length(dNodesFun)>0
    dBCarr = [dNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in dNodes]
    push!(bcValNodes,dBCarr)
    push!(bcValNamesSym,:dBC)
  end
  dNodesArr::Vector{DirichletArray} = filter(isDirichletArray,bcfun)
  if length(dNodesArr)>0
    dBCarr = [dNodesArr[1].v[i] for i in dNodes]
    push!(bcValNodes,dBCarr)
    push!(bcValNamesSym,:dBC)
  end

  # neumann
  nNodesFun::Vector{NeumannFunction} = filter(isNeumannFunction,bcfun)
  nNodesArr::Vector{NeumannArray} = filter(isNeumannArray,bcfun)
  if (length(nNodesFun)>0) || (length(nNodesArr)>0)
    if length(nNodesFun)>0
      nBCarr = [nNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                                for i=1:length(nNodes)]
      push!(bcValNodes,nBCarr)
      push!(bcValNamesSym,:nBC)
    elseif length(nNodesArr)>0
      nBCarr = [nNodesArr[1].v[i] for i=1:length(nNodes)]
      push!(bcValNodes,nBCarr)
      push!(bcValNamesSym,:nBC)
    end
  else
    if :nNodes in bcIDNamesSym
      nBCarr = zeros(Float64,length(nNodes))
      push!(bcValNodes,nBCarr)
      push!(bcValNamesSym,:nBC)
    end
  end

  # robin
  rNodesFun::Vector{RobinFunction} = filter(isRobinFunction,bcfun)
  if length(rNodesFun)>0
    rBCarr = [rNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in rNodes]
    push!(bcValNodes,rBCarr)
    push!(bcValNamesSym,:rBC)
  end
  rNodesArr::Vector{DirichletArray} = filter(isRobinArray,bcfun)
  if length(rNodesArr)>0
    rBCarr = [rNodesArr[1].v[i] for i in rNodes]
    push!(bcValNodes,rBCarr)
    push!(bcValNamesSym,:rBC)
  end

  # forcing
  fNodesFun::Vector{ForcingFunction} = filter(isForcingFunction,bcfun)
  fNodesArr::Vector{ForcingArray} = filter(isForcingArray,bcfun)
  if (length(fNodesFun)>0) || (length(fNodesArr)>0)
    if length(fNodesFun)>0
      farr = [fNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                                for i=1:length(mesh.xy)]
      push!(bcValNodes,farr)
      push!(bcValNamesSym,:forcing)
    elseif length(fNodesArr)>0
      farr = [fNodesArr[1].v[i] for i=1:length(mesh.xy)]
      push!(bcValNodes,farr)
      push!(bcValNamesSym,:forcing)
    end
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
  return Problem(OpType,bcid,bcval)
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
    dUBCarr = [dNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in dUNodes]
    push!(bcValNodes,dUBCarr)
    push!(bcValNamesSym,:dUBC)

    dVBCarr = [dNodesFun[2].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i in dVNodes]
    push!(bcValNodes,dVBCarr)
    push!(bcValNamesSym,:dVBC)
  end
  dNodesArr::Vector{DirichletFunction} = filter(isDirichletArray,bcfun)
  if length(dNodesArr)==2
    dUBCarr = [dNodesArr[1].v[i] for i in dUNodes]
    push!(bcValNodes,dUBCarr)
    push!(bcValNamesSym,:dUBC)

    dVBCarr = [dNodesArr[2].v[i] for i in dVNodes]
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
  fNodesArr::Vector{ForcingArray} = filter(isForcingArray,bcfun)
  if length(fNodesFun)==2
    farrX = [fNodesFun[1].f(mesh.xy[i].x,mesh.xy[i].y)
                               for i=1:length(mesh.xy)]
    farrY = [fNodesFun[2].f(node.x,node.y)
                               for node in mesh.xy]
    push!(bcValNodes,farrX)
    push!(bcValNamesSym,:forcingX)
    push!(bcValNodes,farrY)
    push!(bcValNamesSym,:forcingY)
  elseif length(fNodesArr)==2
    farrX = [fNodesArr[1].v[i] for i=1:length(mesh.xy)]
    farrY = [fNodesArr[2].v[i] for i=1:length(mesh.xy)]
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

  # return Problem object
  return Problem(OperatorType,bcid,bcval)
end

####### Problem Methods #######

bc(prob::Problem) = indexlist(prob.bcid)
