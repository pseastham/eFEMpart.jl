module eFEMpart

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
       Mesh

# MeshTransform.jl
export pointTransform,
       onSegment,
       doIntersect,
       isInsideDomain,
       NodesToArray,
       ElementsToArray,
       meshTransform

# Basis.jl
export derivShape2D,
       shapeEval

# MatrixGeneration.jl
export Mass2DMatrix

# TimeStepping.jl
export FEMForwardEuler

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

# ======================================================
# FILES TO LOAD
# ======================================================

# Custom types
include("MeshTypes.jl")
include("ParameterTypes.jl")
include("ProblemTypes.jl")
include("SolutionTypes.jl")

# Mesh input/generation
include("MeshGeneration.jl")
include("GMSHreader.jl")
include("MeshTransform.jl")

# Generate Matrices
include("Basis.jl")
include("BilinearForms.jl")
include("MatrixGeneration.jl")
include("BoundaryConditions.jl")

# Solve Problem
include("Assembly.jl")
include("TimeStepping.jl")
include("Solver.jl")

# PostProcessing
include("vtkExport.jl")
include("ErrorAnalysis.jl")
include("DomainQuadrature.jl")
include("PostProcessing.jl")
include("TracerGenerate.jl")

end # module eFEMpart
