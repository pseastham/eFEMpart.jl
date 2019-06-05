# solver file to be loaded into eFEMpart

using IterativeSolvers, Preconditioners, LinearMaps

################################
###### General Framework #######
################################

# General solver for scalar problems
function solve(prob::Problem,mesh::S,param::T) where 
                {S<:AbstractMesh,T<:AbstractParameter}
  # generate linear operator matrix
  LinOp = GenerateSystem(mesh,prob,param)

  # Apply boundary conditions
  ApplyBC!(LinOp,mesh,prob,param,prob.OperatorType)

  # apply solver
  sol = LinearSolve(mesh,LinOp,prob,param,prob.OperatorType)

  return sol
end
# using new boundary condition application during matrix assembly
function solve_darcyvar_2b(prob,mesh,param)
    LinOp = LinearOperator(mesh,prob)

    # generate linear operator matrix & apply BC!
    start = time()
    LinOp.Op,LinOp.rhs = Darcy2DMatrix_2b(mesh,prob.bcval[:forcing],param,prob)
    #println("linear system generation & bc apply (together):")
    println("linear system generation & bc apply:")
    println("  ",time() - start," seconds")

    #start = time()
    #ApplyBC!(LinOp,mesh,prob,param,prob.OperatorType)
    #println("linear system generation & bc apply (together):")
    #println("  ",time() - start," seconds")

    # apply solver
    start = time()
    #stemp = @elapsed(sol = LinearSolve(mesh,LinOp,prob,param,prob.OperatorType))
    sol = LinearSolve(mesh,LinOp,prob,param,prob.OperatorType)
    println("solver time:")
    println("  ",time() - start," seconds")

    return sol
end

# solve function without parameter included (for :Poisson2D)
function solve(prob::Problem,mesh::S) where S<:AbstractMesh
   param = PoissonParam(1.0);  return solve(prob,mesh,param)
end

# General solver for scalar problems
function solve!(sol::R,prob::Problem,mesh::S,param::T) where 
                {R<:AbstractSolution,S<:AbstractMesh,T<:AbstractParameter}
  # generate linear operator matrix
  LinOp = GenerateSystem(mesh,prob,param)

  # Apply boundary conditions
  ApplyBC!(LinOp,mesh,prob,param,prob.OperatorType)

  # apply solver
  LinearSolve!(sol,mesh,LinOp,prob,param,prob.OperatorType)

  return nothing
end

##############################
###### Generate Matrix #######
##############################

function GenerateSystem(mesh::R,prob::S,param::T) where 
          {R<:AbstractMesh,S<:AbstractProblem,T<:AbstractParameter}
  LinOp = LinearOperator(mesh,prob)

  # POISSON ✓
  if prob.OperatorType == :Poisson2D
    LinOp.Op,LinOp.rhs = Laplace2DMatrix(mesh,prob.bcval[:forcing],param)

  # DARCY -- NOT TESTED YET
  elseif prob.OperatorType == :Darcy2D
    LinOp.Op,LinOp.rhs = Darcy2DMatrix(mesh,prob.bcval[:forcing],param)

  # ADVECTION-DIFFUSION 2D ✓
  elseif prob.OperatorType == :AdvDiff2D
    LinOp.Op,LinOp.rhs = AdvDiff2DMatrix(mesh,prob.bcval[:forcing],param)

  # Stokes 2D ✓
  elseif prob.OperatorType == :Stokes2D
    LinOp.Op,LinOp.rhs = Stokes2DMatrix(mesh,prob,param)

  # BRINKMAN 2D  -- NOT TESTED YET
  elseif prob.OperatorType == :Brinkman2D
    LinOp.Op,LinOp.rhs = Brinkman2DMatrix(mesh,prob,param)

  # BRINKMAN MULTIPHASE 2D -- NOT TESTED YET
  elseif prob.OperatorType == :BrinkmanMP2D
    LinOp.Op,LinOp.rhs = BrinkmanMP2DMatrix(mesh,prob,param)

  # ADVECTION-DIFFUSION AXISYMMETRIC
  elseif prob.OperatorType == :AdvDiffAS
    LinOp.Op,LinOp.rhs = AdvDiffASMatrix(mesh,prob.bcval[:forcing],param)

  # STOKES AXISYMMETRIC
  elseif prob.OperatorType == :StokesAS
    LinOp.Op,LinOp.rhs = StokesASMatrix(mesh,prob,param)

  end

  return LinOp
end

########################################
###### Apply Boundary Conditions #######
########################################

# solve boundary condition
function ApplyBC!(LinOp::LinearOperator,mesh::S,prob::Problem,param::T,
                  OpType::Symbol) where {S<:AbstractMesh,T<:AbstractParameter}
  # POISSON 2D or DARCY 2D
  if OpType == :Poisson2D || OpType == :Darcy2D
    if :dNodes in bc(prob)
      scalarDirichlet!(prob.bcid[:dNodes],prob.bcval[:dBC],LinOp.Op,LinOp.rhs)
    end
    if :nNodes in bc(prob)
      scalarNeumann!(mesh.xy,mesh.cm,mesh.order,
                     prob.bcid[:nNodes],prob.bcval[:nBC],LinOp.rhs)
    end

  # ADVECTION-DIFFUSION 2D
  elseif OpType == :AdvDiff2D
    if :dNodes in bc(prob)
      scalarDirichlet!(prob.bcid[:dNodes],prob.bcval[:dBC],LinOp.Op,LinOp.rhs)
    end

    if :nNodes in bc(prob)
      scalarNeumann!(mesh.xy,mesh.cm,mesh.order,
                     prob.bcid[:nNodes],prob.bcval[:nBC],LinOp.rhs)
    end

    if :rNodes in bc(prob)
      AdvDiffRobin!(mesh.xy,mesh.cm,prob.bcid[:rNodes],
                    LinOp.Op,LinOp.rhs,hcat(param.wx,param.wy))
    end

  # STOKES 2D or BRINKMAN 2D or BRINKMAN MULTIPHASE 2D
  elseif OpType == :Stokes2D || OpType == :Brinkman2D || OpType == :BrinkmanMP2D
    if :dUNodes in bc(prob)
      fluidDirichlet!(mesh,prob,LinOp)
    end

    if !(:nNodes in bc(prob))
      LinOp.Op,LinOp.rhs = fluidNeumann(mesh,prob,LinOp)
    end

  # ADVECTION-DIFFUSION AXISYMMETRIC
  elseif OpType == :AdvDiffAS
    if :dNodes in bc(prob)
      scalarDirichlet!(prob.bcid[:dNodes],prob.bcval[:dBC],LinOp.Op,LinOp.rhs)
    end

    if :nNodes in bc(prob)
      scalarNeumann!(mesh.xy,mesh.cm,mesh.order,
                    prob.bcid[:nNodes],prob.bcval[:nBC],LinOp.rhs)
    end

    if :rNodes in bc(prob)
      AdvDiffRobin!(mesh.xy,mesh.cm,prob.bcid[:rNodes],
                    LinOp.Op,LinOp.rhs,prob.wind)
    end

  # STOKES AXISYMMETRIC  
  elseif OpType == :StokesAS
    if :dUNodes in bc(prob)
      fluidDirichlet!(mesh,prob,LinOp)
    end

    if !(:nNodes in bc(prob))
      LinOp.Op,LinOp.rhs = fluidNeumann(mesh,prob,LinOp)
    end
  end
end

##############################
######## Solve System ########
##############################

function LinearSolve!(sol::AbstractSolution,mesh::AbstractMesh,
                      LinOp::LinearOperator,prob::AbstractProblem,
                      param::AbstractParameter,OpType::Symbol)
  # POISSON 2D                  
  if OpType == :Poisson2D
    sol.u = GaussElimSolve(LinOp)

  # DARCY 2D 
  elseif OpType == :Darcy2D
    println("LinearSolve! not defined for Darcy")
    return FluidSolution(u,v,p)

  # ADVECTION-DIFFUSION 2D   
  elseif OpType == :AdvDiff2D || OpType == :AdvDiffAS
    sol.u = GaussElimSolve(LinOp)
  
  # STOKES 2D or BRINKMAN 2D or BRINKMAN MULTIPHASE 2D
  elseif ( OpType == :Stokes2D     || OpType == :Brinkman2D || 
           OpType == :BrinkmanMP2D)
    fluidSolve!(sol,mesh,prob,LinOp)

  # STOKES AS requires row-normalization preconditioner!
  elseif OpType == :StokesAS
    #rowNormalizePreconditioner!(LinOp)
    fluidSolve!(sol,mesh,prob,LinOp)
  end

  return nothing
end

function LinearSolve(mesh::AbstractMesh,
                      LinOp::LinearOperator,prob::AbstractProblem,
                      param::AbstractParameter,OpType::Symbol)
  # POISSON 2D                  
  if OpType == :Poisson2D
    u = GaussElimSolve(LinOp)
    #u = cg(LinOp.Op,LinOp.rhs)

  return ScalarSolution(u)

  # DARCY 2D 
  elseif OpType == :Darcy2D
    # preconditioner=
    pc = AMGPreconditioner{SmoothedAggregation}(LinOp.Op)

    p = IterativeSolvers.cg(LinOp.Op, LinOp.rhs;verbose=false,Pl=pc)

    #p = GaussElimSolve(LinOp)
    u,v = DarcyVelocity(mesh,param.α,p)
  return FluidSolution(u,v,p)

  # ADVECTION-DIFFUSION 2D   
  elseif OpType == :AdvDiff2D || OpType == :AdvDiffAS
    #pc = DiagonalPreconditioner(LinOp.Op)

    #u = gmres(LinOp.Op, LinOp.rhs;verbose=true,Pl=pc)
    u = GaussElimSolve(LinOp)
  return ScalarSolution(u)

  # STOKES 2D or BRINKMAN 2D or BRINKMAN MULTIPHASE 2D
  elseif ( OpType == :Stokes2D     || OpType == :Brinkman2D || 
    OpType == :BrinkmanMP2D)
    return fluidSolve(mesh,prob,LinOp)
  
  # STOKES AS requires row-normalization preconditioner!
  elseif OpType == :StokesAS
    #rowNormalizePreconditioner!(LinOp)
    return fluidSolve(mesh,prob,LinOp)
  end
end

function GaussElimSolve(LinOp::LinearOperator)
  U = LinOp.Op\LinOp.rhs

  return U
end

"""
  fluidSolve(xy,cm,nNodes,Stiff,F;PRINT=false) -> u,v,p

Computes and decomposes the array-solution U into vector velocity [u,v] and scalar pressure p
"""
function fluidSolve(mesh::FluidMesh,prob::Problem,
                    LinOp::AbstractLinearOperator;PRINT=false)
  # preconditioner
  #pc = IterativeSolvers.Identity()
  #pc   = Preconditioners.DiagonalPreconditioner(LinOp.Op)
  #soln = IterativeSolvers.minres(LinOp.Op, LinOp.rhs;
  #                              verbose=true)

  soln = qr(LinOp.Op)\(LinOp.rhs)
  #soln = lu(LinOp.Op)\(LinOp.rhs)

  # decompose solution
  NumUnodes = length(mesh.xy)
  u = soln[1:NumUnodes]
  v = soln[NumUnodes+1:2*NumUnodes]
  p = soln[2*NumUnodes+1:end]

  # adjust for enclosed boundary condition
  if !(:nNodes in bc(prob))
    if PRINT
      println("adjusting for enclosed boundary condition")
    end
    p=p[1:end-1]
  end

  return FluidSolution(u,v,p)
end

"""
  fluidSolve(xy,cm,nNodes,Stiff,F;PRINT=false) -> u,v,p

Computes and decomposes the array-solution U into vector velocity [u,v] and scalar pressure p
"""
function fluidSolve!(sol::FluidSolution,
                     mesh::FluidMesh,prob::Problem,
                     LinOp::AbstractLinearOperator;PRINT=false)
  soln = lu(LinOp.Op)\LinOp.rhs

  # decompose solution
  NumUnodes = length(mesh.xy)
  sol.u = soln[1:NumUnodes]
  sol.v = soln[NumUnodes+1:2*NumUnodes]
  sol.p = soln[2*NumUnodes+1:end]

  # adjust for enclosed boundary condition
  if !(:nNodes in bc(prob))
    if PRINT
      println("adjusting for enclosed boundary condition")
    end
    sol.p=sol.p[1:end-1]
  end

  return nothing
end

function DarcyVelocity(mesh,alpha::Vector{Float64},p::Vector{Float64})
  Nel    = length(mesh.cm)
  Nnodes = length(mesh.xy)

  nodeCounter = zeros(Int,Nnodes)

  u = zeros(Float64,Nnodes)
  v = zeros(Float64,Nnodes)

  xNodes = zeros(Float64,9)
  yNodes = zeros(Float64,9)
  αNodes = zeros(Float64,9)
  pNodes = zeros(Float64,9)

  w,s,t = GaussQuadPoints2D(3)

  # for each element
  for el=1:Nel
    c = mesh.cm[el].NodeList
    for i=1:9
      xNodes[i] = mesh.xy[c[i]].x
      yNodes[i] = mesh.xy[c[i]].y
      αNodes[i] = alpha[c[i]]
      pNodes[i] = p[c[i]]
    end

    # for each node in element (assumes quadratic basis)
    for i=1:9
      phi,dphidx,dphidy,jac = derivShape2D(s[i],t[i],xNodes,yNodes,2)
      dpdx = shapeEval(pNodes,dphidx)
      dpdy = shapeEval(pNodes,dphidy)
      αg = shapeEval(αNodes,phi)

      u[c[i]] += -αg*dpdx/jac
      v[c[i]] += -αg*dpdy/jac
    end
    nodeCounter[c] .+= 1
  end

  # take average
  for i=1:Nnodes
    u[i] /= nodeCounter[i]
    v[i] /= nodeCounter[i]
  end

  return u,v
end
function DarcyVelocity(mesh,alpha::Float64,p::Vector{Float64})
  return DarcyVelocity(mesh,alpha*ones(Float64,length(mesh.xy)),p)
end

function rowNormalizePreconditioner!(LinOp)
    # apply preconditioner
    rows,cols,_ = findnz(LinOp.Op)
    nrows = length(rows)

    # compute row norms
    counter = 0
    alreadyDone = []
    for ti=1:nrows
      if rows[ti] in alreadyDone

      else
        normval = norm(LinOp.Op[rows[ti],:])
        
        # normalize Matrix
        LinOp.Op[rows[ti],:] ./= normval
    
        # normalize forcing -- this is going twice!!!
        LinOp.rhs[rows[ti]] /= normval

        push!(alreadyDone,rows[ti])
        counter += 1
      end
    end

  return nothing
end