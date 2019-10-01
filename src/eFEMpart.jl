module eFEMpart

# for access to SparseMatrixCSC
using SparseArrays
# for access to lufact
using LinearAlgebra
 # particularly used in Solver.jl
using IterativeSolvers, Preconditioners
# used for testing
using BenchmarkTools

# ======================================================
# FILES TO LOAD
# ======================================================

# Custom types
include("eFEM/MeshTypes.jl")
include("eFEM/ParameterTypes.jl")
include("eFEM/ProblemTypes.jl")
include("eFEM/SolutionTypes.jl")

# Mesh input/generation
include("eFEM/MeshGeneration.jl")
include("eFEM/GMSHreader.jl")
include("eFEM/MeshTransform.jl")

# Generate Matrices
include("eFEM/Basis.jl")
include("eFEM/BilinearForms.jl")
include("eFEM/MatrixGeneration.jl")
include("eFEM/BoundaryConditions.jl")

# Solve Problem
include("eFEM/Assembly.jl")
include("eFEM/TimeStepping.jl")
include("eFEM/Solver.jl")

# PostProcessing
include("eFEM/vtkExport.jl")
include("eFEM/ErrorAnalysis.jl")
include("eFEM/DomainQuadrature.jl")
include("eFEM/PostProcessing.jl")
include("eFEM/TracerGenerate.jl")

# Particles
include("part/CellLists.jl")
include("part/fgt.jl")
include("part/ParticleTypes.jl")
include("part/UpdateParticles.jl")

include("part/AdhesionForce.jl")

# ======================================================
# EXPORTED FUNCTIONS
# ======================================================

# MeshTypes.jl
export AbstractMesh

# ParameterTypes.jl
export PoissonParam,
       HeatParam,
       AdvDiffParam,
       DarcyParam,
       FluidParam,
       BrinkmanParam,
       BrinkmanMPParam,
       FluidVarParam,
       AbstractVariableParameter,
       AbstractConstantParameter

# ProblemTypes.jl
export Problem,
       Dirichlet,
       Neumann,
       Robin,
       Forcing

# SolutionsTypes.jl
export ScalarSolution,
       FluidSolution

# MeshGeneration.jl
export squareMesh,
       squareMeshFluid,
       axisymMesh,
       axisymFluid,
       FluidMesh,
       ScalarMesh,
       Mesh,
       meshgrid,
       ElementsToArray,
       NodesToArray

# MeshTransform.jl
export pointTransform,
       onSegment,
       doIntersect,
       isInsideDomain,
       meshTransform

# Basis.jl
export derivShape2D,
       shapeEval,
       GaussEdge

# MatrixGeneration.jl
export Mass2DMatrix

# TimeStepping.jl
export FEMForwardEuler,
       progressBar,
       sToHMS

# Solver.jl
export solve,
       solve!,
       GenerateSystem,
       ApplyBC!

# vtkExport.jl
export vtksave,
       Path,
       ScalarData,
       ScalarNames,
       VectorData,
       VectorNames

# ErrorAnalysis.jl
export DomainNorm,
       hCalc

# DomainQuadrature.jl
export DomainQuad,
       SurfaceQuad,
       SurfaceFlux,
       Surface,
       VelocityFlux

# PostProcessing.jl
export computeStress,
       computeDerivative,
       computeExtension,
       computeCompressibility

# TracerGenerate.jl
export GenerateTracers,
       Tracer,
       TracerInfo,
       TracerLineSource

# -----
# part
# -----

export AbstractWall, NearestPoint, generateQuadNodes!,isInLine

# CellLists.jl
export femGenerateMap,
       generateCellList,
       updateCellList!

# fgt.jl
export interpFGT!

# ParticleTypes.jl
export Point2D,
       LineWall,
       CircleWall,
       ArcWall

# UpdateParticles.jl
export updateParticle_all!,
       updateParticle_all_nofluid!

# for testing...
export BarycentricVelocityInterp_CL!,
       computeCohesion!,
       AdhesionForce!,
       computeCohesion_CL!,
       LennardJonesForceMagnitude

end # module eFEMpart
