module eFEMpart

# for access to SparseMatrixCSC
using SparseArrays
# for access to lufact
using LinearAlgebra

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
       shapeEval

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

# SinkholeTypes.jl
export AbstractWall,            # Abstract type that encapsulates all concrete wall types
       LineWall,                # mutable type for describing line walls
       CircleWall,              # mutable type for describing circle walls
       ArcWall,                 # mutable type for describing circle arc walls
       Particle,
       Point,                   # mutable type for points
       ParticleList             # Initializes list of particles

# CellLists.jl
export Cell,
       CellGenerate      # Generates cell lists from either a list
                         #   of particles or a mesh

# Porosity.jl
export SandToPorosityGaussian

# ParticleForces.jl
export LennardJonesForce,
       LennardJonesPotential,
       LennardJonesPotentialMagnitude

# WallRepulsion.jl 
export NearestPoint,
       isCloseEnough,
       GenerateQuadNodes,
       isInLine,
       WallTrapQuad

# UpdateParticles.jl
export UpdateParticle_all!,
       UpdateParticle_novelocity!

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
include("part/SinkholeTypes.jl")
include("part/Walls.jl")
include("part/CellLists.jl")
include("part/Porosity.jl")
include("part/ParticleForces.jl")
include("part/WallRepulsion.jl")
include("part/UpdateParticles.jl")

end # module eFEMpart
