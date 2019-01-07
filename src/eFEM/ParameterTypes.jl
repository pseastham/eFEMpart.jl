# File loaded into eFEM

### Abstract Types ###

abstract type AbstractParameter end

abstract type AbstractDummyParameter <: AbstractParameter end
abstract type AbstractConstantParameter <: AbstractParameter end
abstract type AbstractVariableParameter <: AbstractParameter end

### Concrete Types Definitions ###

mutable struct PoissonParam <: AbstractDummyParameter
  κ::Float64
end

mutable struct HeatConstParam <: AbstractConstantParameter
  κ::Float64
end
mutable struct HeatVarParam <: AbstractVariableParameter
  κ::Vector{Float64}
end

mutable struct DarcyConstParam <: AbstractConstantParameter
  α::Float64
end
mutable struct DarcyVarParam <: AbstractVariableParameter
  α::Vector{Float64}
end

mutable struct AdvDiffConstParam <: AbstractConstantParameter
  u::Vector{Float64}
  v::Vector{Float64}
  κ::Float64
end
mutable struct AdvDiffVarParam <: AbstractVariableParameter
  u::Vector{Float64}
  v::Vector{Float64}
  κ::Vector{Float64}
end

mutable struct FluidConstParam <: AbstractConstantParameter
  μ::Float64
end
mutable struct FluidVarParam <: AbstractVariableParameter
  μ::Vector{Float64}
end

mutable struct BrinkmanConstParam <: AbstractConstantParameter
  μ::Float64
  α::Float64
end
mutable struct BrinkmanVarParam <: AbstractVariableParameter
  μ::Vector{Float64}
  α::Vector{Float64}
end

mutable struct BrinkmanMPConstParam <: AbstractConstantParameter
  α1::Float64
  α2::Float64
  α3::Float64
end
mutable struct BrinkmanMPVarParam <: AbstractVariableParameter
  α1::Vector{Float64}
  α2::Vector{Float64}
  α3::Vector{Float64}
end

### Abstract Problem-specific Unions ###

AbstractHeat       = Union{HeatConstParam,HeatVarParam}
AbstractAdvDiff    = Union{AdvDiffConstParam,AdvDiffVarParam}
AbstractDarcy      = Union{DarcyConstParam,DarcyVarParam}
AbstractFluid      = Union{FluidConstParam,FluidVarParam}
AbstractBrinkman   = Union{BrinkmanConstParam,BrinkmanVarParam}
AbstractBrinkmanMP = Union{BrinkmanMPConstParam,BrinkmanMPVarParam}

### Constructor Definitions ###

function HeatParam(κ::T) where T<:Real
  return HeatConstParam(κ)
end
function HeatParam(κ::Vector{T}) where T<:Real
  return HeatVarParam(κ)
end

function AdvDiffParam(u,v,κ::T) where T<:Real
  return AdvDiffConstParam(u,v,κ)
end
function AdvDiffParam(u,v,κ::Vector{T}) where T<:Real
  return AdvDiffVarParam(u,v,κ)
end

function DarcyParam(α::T) where T<:Real
  return DarcyConstParam(α)
end
function DarcyParam(α::Vector{T}) where T<:Real
  return DarcyVarParam(α)
end

function FluidParam(μ::T) where T<:Real
  return FluidConstParam(μ)
end
function FluidParam(μ::Vector{T}) where T<:Real
  return FluidVarParam(μ)
end

function BrinkmanParam(μ::T,α::T) where T<:Real
  return BrinkmanConstParam(μ,α)
end
function BrinkmanParam(μ::T,α::Vector{T}) where T<:Real
  return BrinkmanVarParam(μ*zeros(T,length(α)),α)
end
function BrinkmanParam(μ::Vector{T},α::T) where T<:Real
  return BrinkmanVarParam(μ,α*ones(T,length(μ)))
end
function BrinkmanParam(μ::Vector{T},α::Vector{T}) where T<:Real
  return BrinkmanVarParam(μ,α)
end

function BrinkmanMPParam(α1::T,α2::T,α3::T) where T<:Real
  return BrinkmanMPConstParam(α1,α2,α3)
end
function BrinkmanMPParam(α1::Vector{T},α2::T,α3::T) where T<:Real
  o = ones(T,length(α1))
  return BrinkmanMPVarParam(α1,α2*o,α3*o)
end
function BrinkmanMPParam(α1::T,α2::Vector{T},α3::T) where T<:Real
  o = ones(T,length(α2))
  return BrinkmanMPVarParam(α1*o,α2,α3*o)
end
function BrinkmanMPParam(α1::T,α2::T,α3::Vector{T}) where T<:Real
  o = ones(T,length(α1))
  return BrinkmanMPVarParam(α1,α2*o,α3*o)
end
function BrinkmanMPParam(α1::Vector{T},α2::Vector{T},α3::T) where T<:Real
  o = ones(T,length(α1))
  return BrinkmanMPVarParam(α1,α2,α3*o)
end
function BrinkmanMPParam(α1::Vector{T},α2::T,α3::Vector{T}) where T<:Real
  o = ones(T,length(α1))
  return BrinkmanMPVarParam(α1,α2*o,α3)
end
function BrinkmanMPParam(α1::T,α2::Vector{T},α3::Vector{T}) where T<:Real
  o = ones(T,length(α2))
  return BrinkmanMPVarParam(α1*o,α2,α3)
end
function BrinkmanMPParam(α1::Vector{T},α2::Vector{T},α3::Vector{T}) where T<:Real
  return BrinkmanMPVarParam(α1,α2,α3)
end