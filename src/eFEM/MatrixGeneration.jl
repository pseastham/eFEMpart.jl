# functions to be loaded into eFEMpart

####################################
###### Application Functions #######
####################################

# MASS
function Mass2DMatrix(mesh)
    return assembleScalar(mesh,localMass2D!,0.0)
end

# LAPLACE 2D
function Laplace2DMatrix(mesh,farr,param::T) where T<:AbstractDummyParameter
    D = assembleScalar(mesh,localLaplace2D!,0.0)
    F = WeakScalar2D(mesh,farr)
    return D,F 
end

# DARCY 2D
function Darcy2DMatrix(mesh,farr,param::T) where T<:AbstractConstantParameter
    D = assembleScalar(mesh,localLaplace2D!,0.0)
    D = param.α*D
    F = WeakScalar2D(mesh,farr)
    return D,F 
end
function Darcy2DMatrix(mesh,farr,param::T) where T<:AbstractVariableParameter
    D = assembleScalar(mesh,localDarcy2D!,param.α)
    F = WeakScalar2D(mesh,farr)
    return D,F 
end

# ADVECTION-DIFFUSION 2D
function AdvDiff2DMatrix(mesh,farr,param::T) where T<:AbstractConstantParameter
    D = assembleScalar(mesh,localLaplace2D!,0.0)
    A = assembleScalar(mesh,localAdvDiff2D!,param)
    Stiff = A + param.κ*D
    F = WeakScalar2D(mesh,farr)
    return Stiff, F
end
function AdvDiff2DMatrix(mesh,farr,param::T) where T<:AbstractVariableParameter
    D = assembleScalar(mesh,localDarcy2D!,param.κ)
    A = assembleScalar(mesh,localAdvDiff2D!,param)
    Stiff = A + D
    F = WeakScalar2D(mesh,farr)
    return Stiff, F
end

# STOKES 2D
function Stokes2DMatrix(mesh::FluidMesh,prob,param::T) where T<:AbstractConstantParameter
    Stiff = assembleHalfFluid(mesh,localStokesConst2D!,param.μ)
    Fx    = WeakScalar2D(mesh,prob.bcval[:forcingX])
    Fy    = WeakScalar2D(mesh,prob.bcval[:forcingY])
    G     = zeros(length(mesh.xyp))
    F     = vcat(Fx,Fy,G)

    return Stiff, F
end
function Stokes2DMatrix(mesh::FluidMesh,prob,param::T) where T<:AbstractVariableParameter
    # new original -- not working
    #Stiff = assembleFullFluid(mesh,localStokesVar2D!,param.μ)
    
    # new alternative -- not working
    #Stiff = assembleFullFluid(mesh,localStokesAlternateVar2D!,param.μ)
    
    # old original -- THIS ONE IS WORKING
    Stiff = assembleFullFluid_WRONG(mesh,localStokesVar2D_WRONG!,param.μ)
    
    Fx    = WeakScalar2D(mesh,prob.bcval[:forcingX])
    Fy    = WeakScalar2D(mesh,prob.bcval[:forcingY])
    G     = zeros(length(mesh.xyp))
    F     = vcat(Fx,Fy,G)

    return Stiff, F
end

# BRINKMAN 2D
function Brinkman2DMatrix(mesh::FluidMesh,prob,param::T) where T<:AbstractConstantParameter
    Stiff = assembleHalfFluid(mesh,localBrinkmanConst2D!,param.μ,param.α)
    Fx    = WeakScalar2D(mesh,prob.bcval[:forcingX])
    Fy    = WeakScalar2D(mesh,prob.bcval[:forcingY])
    G     = zeros(length(mesh.xyp))
    F     = vcat(Fx,Fy,G)

    return Stiff, F
end
function Brinkman2DMatrix(mesh::FluidMesh,prob,param::T) where T<:AbstractVariableParameter
    Stiff = assembleFullFluid(mesh,localBrinkmanVar2D!,param.μ,param.α)
    Fx    = WeakScalar2D(mesh,prob.bcval[:forcingX])
    Fy    = WeakScalar2D(mesh,prob.bcval[:forcingY])
    G     = zeros(length(mesh.xyp))
    F     = vcat(Fx,Fy,G)

    return Stiff, F
end

"""
Generates matrix operator to be used to solve brinkman-like equations
  -α1*lap(u) + α2*u + α3*grad(p) = f
"""
function BrinkmanMP2DMatrix(mesh::FluidMesh,prob,param::T) where T<:AbstractConstantParameter
	Stiff = assembleHalfFluid(mesh,localBrinkmanMPConst2D!,param.α1,param.α2,param.α3)
	Fx    = WeakScalar2D(mesh,prob.bcval[:forcingX])
	Fy    = WeakScalar2D(mesh,prob.bcval[:forcingY])
	G     = zeros(length(mesh.xyp))
	F     = vcat(Fx,Fy,G)

	return Stiff, F
end
function BrinkmanMP2DMatrix(mesh::FluidMesh,prob,param::T) where T<:AbstractVariableParameter
    #Stiff = assembleFullFluid(mesh,localMPBVar2D!,param.α1,param.α2,param.α3)
    Stiff = assembleMPB(mesh,localMPBVar2D!,param.α1,param.α2,param.α3)
	Fx    = WeakScalar2D(mesh,prob.bcval[:forcingX])
	Fy    = WeakScalar2D(mesh,prob.bcval[:forcingY])
	G     = zeros(length(mesh.xyp))
	F     = vcat(Fx,Fy,G)

	return Stiff, F
end

# ADVECTION-DIFFUSION AXISYMMETRIC
function AdvDiffASMatrix(mesh,farr,param::T) where T<:AbstractConstantParameter
    D = assembleScalar(mesh,localLaplaceConstAS!,0.0)
    A = assembleScalar(mesh,localAdvDiffAS!,param)
    Stiff = A + param.κ*D
    F = WeakScalarAS(mesh,farr)
    return Stiff, F
end
function AdvDiffASMatrix(mesh,farr,param::T) where T<:AbstractVariableParameter
    D = assembleScalar(mesh,localLaplaceVarAS!,param.κ)
    A = assembleScalar(mesh,localAdvDiffAS!,param)
    Stiff = A + D
    F = WeakScalarAS(mesh,farr)
    return Stiff, F
end

# STOKES AXISYMMETRIC
function StokesASMatrix(mesh,prob,param::T) where T<:AbstractConstantParameter
    Stiff = assembleHalfFluid(mesh,localStokesConstAS!,param.μ)
    Fx    = WeakScalarAS(mesh,prob.bcval[:forcingX])
    Fy    = WeakScalarAS(mesh,prob.bcval[:forcingY])
    G     = zeros(Float64,length(mesh.xyp))
    F     = vcat(Fx,Fy,G)

    return Stiff, F   
end
function StokesASMatrix(mesh,prob,param::T) where T<:AbstractVariableParameter
    Stiff = assembleFullFluid_WRONG(mesh,localStokesVarAS!,param.μ)
    #Stiff = assembleHalfFluid(mesh,localStokesVarASAlternate!,param.μ)
    Fx    = WeakScalarAS(mesh,prob.bcval[:forcingX])
    Fy    = WeakScalarAS(mesh,prob.bcval[:forcingY])
    G     = zeros(Float64,length(mesh.xyp))
    F     = vcat(Fx,Fy,G)

    return Stiff, F   
end