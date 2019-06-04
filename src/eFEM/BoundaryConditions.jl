# File loaded into eFEMpart

"""
  scalarDirichlet!(dirichletNodes,dBCarr,Stiff,F)

Applies Dirichlet boundary conditions for Scalar equation
"""
function scalarDirichlet!(dirichletNodes,dBCarr,Stiff,F)
  Nnodes  = length(F)
  dLength = length(dirichletNodes)
  rvals = Stiff.rowval

  for i=1:dLength
    d = dirichletNodes[i]

    for j=1:Nnodes
    #for j in rvals
      # adjust RHS with dirichlet BC
	    F[j] -= Stiff[j,d]*dBCarr[i]

	    # zero out (dirichlet) BC locations
	    Stiff[j,d] = 0.0
	    Stiff[d,j] = 0.0
    end

    #fill!(view(Stiff,:,d),0.0)
    #fill!(view(Stiff,d,:),0.0)

    # identity in correct location
    Stiff[d,d] = 1.0

    # set dirichlet nodes to exact boundary condition
    F[d] = dBCarr[i]
  end
end


"""
  scalarNeumann!(xy,cm,order,neumannNodes,nBCarr,F)

Applies Neumann boundary conditions for arbitrary scalar equation
"""
function scalarNeumann!(xy,cm,order,neumannNodes,nBCarr,F)
  Nel = size(cm,1)
  gaussorder=0

  # generate basis values at gauss points
  if order == :Linear
    order=1
    gaussorder = 2
  elseif order == :Quadratic
    order=2
    gaussorder = 3
  end

  if length(neumannNodes)>0
    for el=1:Nel
      xNodes = [xy[cm[el].NodeList[i]].x for i in 1:Int(gaussorder^2)]
      yNodes = [xy[cm[el].NodeList[i]].y for i in 1:Int(gaussorder^2)]
      # loop over ONLY QUAD NODES of each element
      for i=1:4
        # check that nodes i is a neumann Boundary node
        if cm[el].NodeList[i] in neumannNodes
          # check that next node i+1 is a neumann Boundary node
          j = mod(i,4)+1
          if cm[el].NodeList[j] in neumannNodes
            ci = cm[el].NodeList[i]
            cj = cm[el].NodeList[j]

            w,s,t = GaussEdge(i,gaussorder)

            jac = sqrt((xNodes[i]-xNodes[j])^2+(yNodes[i]-yNodes[j])^2)

            ki = findfirst(s->s==ci,neumannNodes)
            kj = findfirst(s->s==cj,neumannNodes)

            fluxi = nBCarr[ki]
            fluxj = nBCarr[kj]

            for gpt=1:(order+1)
              phi,_,_,_ = derivShape2D(s[gpt],t[gpt],xNodes,yNodes,order)
              fluxVal = shapeEval([fluxi,fluxj],[phi[i],phi[j]])

              F[ci] += w[gpt]*fluxVal*phi[i]*jac
              F[cj] += w[gpt]*fluxVal*phi[j]*jac
            end
          end
        end
      end
    end
  end
end


"""
  scalarBC!(xy,cm,order,dirichletNodes,neumannNodes,
                   dBCarr,nBCarr,Stiff,F)

Applies Dirichlet and Neumann boundary conditions to arbitrary scalar equation
"""
function scalarBC!(xy,cm,order,dirichletNodes,neumannNodes,
                   dBCarr,nBCarr,Stiff,F)
  scalarDirichlet!(dirichletNodes,dBCarr,Stiff,F)
  scalarNeumann!(xy,cm,order,neumannNodes,nBCarr,F)
end


"""
  fluidUDirichlet!(xy,xyp,dUNodes,dUBCarr,Stiff,F)

Applies Dirichlet boundary conditions for horizontal velocity components for Stokes's equation.

Updates the Stiff matrix and F vector
"""
function fluidUDirichlet!(xy,xyp,dUNodes,dUBCarr,Stiff,F)
  nu = length(xy)
  np = length(xyp)
  Nnodes  = length(F)
  uLength = length(dUNodes)

  d = dUNodes
  s1 = 1:nu
  s2 = (nu+1):2*nu
  s3 = (2*nu+1):(2*nu+np)

  s12 = 1:2*nu

  F[s1] -= Stiff[s1,d]*dUBCarr
  F[s2] -= Stiff[s2,d]*dUBCarr
  F[s3] -= Stiff[s3,d]*dUBCarr

  for i=1:uLength
    for j in s12
      Stiff[j,d[i]] = 0.0
      Stiff[d[i],j] = 0.0
    end

    # identity in correct location
    Stiff[d[i],d[i]] = 1.0

    # set dirichlet nodes to exact boundary condition
    F[d[i]] = dUBCarr[i]

    for j in s3
      Stiff[d[i],j] = 0.0
      Stiff[j,d[i]] = 0.0
    end
  end
end


"""
  fluidVDirichlet!(xy,xyp,dVNodes,dVBCarr,Stiff,F)

Applies Dirichlet boundary conditions for vertical velocity components for Stokes's equation.

Updates the Stiff matrix and F vector
"""
function fluidVDirichlet!(xy,xyp,dVNodes,dVBCarr,Stiff,F)
  nu = length(xy)
  np = length(xyp)
  vLength = length(dVNodes)

  d  = nu .+ dVNodes
  s1 = 1:nu
  s2 = (nu+1):2*nu
  s3 = (2*nu+1):(2*nu+np)

  s12 = 1:2*nu

  F[s1] -= Stiff[s1,d]*dVBCarr
  F[s2] -= Stiff[s2,d]*dVBCarr
  F[s3] -= Stiff[s3,d]*dVBCarr

  for i=1:vLength
    for j in s12
      Stiff[j,d[i]] = 0.0
      Stiff[d[i],j] = 0.0
    end

    # identity in correct location
    Stiff[d[i],d[i]] = 1.0

    # set dirichlet nodes to exact boundary condition
    F[d[i]] = dVBCarr[i]

    for j in s3
      Stiff[d[i],j] = 0.0
      Stiff[j,d[i]] = 0.0
    end
  end
end


"""
  fluidDirichlet!(xy,xyp,dUNodes,dVNodes,dUBCarr,dVBCarr,Stiff,F)

computes u and v dirichlet conditions for fluid equation
"""
function fluidDirichlet!(mesh::FluidMesh,prob::Problem,LinOp::LinearOperator)
  # get xy into correct form
  dUNodes = prob.bcid[:dUNodes]
  dUBCarr = prob.bcval[:dUBC]
  dVNodes = prob.bcid[:dVNodes]
  dVBCarr = prob.bcval[:dVBC]

  fluidUDirichlet!(mesh.xy,mesh.xyp,dUNodes,dUBCarr,LinOp.Op,LinOp.rhs)
  fluidVDirichlet!(mesh.xy,mesh.xyp,dVNodes,dVBCarr,LinOp.Op,LinOp.rhs)
end

"""
  fluidNeumann(xy,nNodes,nBCarr,Stiff,F) -> Stiff, F

Because we have assumed nothing else than zero neumann boundary conditions, i.e. 'outflow' conditions, this just checks whether there are no outflow conditions, in which case it adds an additional equation to normalize the pressure about 0.
"""
function fluidNeumann(mesh::FluidMesh,prob::Problem,LinOp::LinearOperator)
  if !(:nNodes in bc(prob))
    NumUnodes = length(mesh.xy)
    NumPnodes = length(mesh.xyp)
    ucol = spzeros(2*NumUnodes,1)
    pcol = ones(NumPnodes,1)/NumPnodes
    urow = ucol'
    prow = pcol'

    s1 = vcat(ucol,pcol)
    s2 = hcat(urow,prow,1)

    LinOp.Op = hcat(LinOp.Op,s1)
    LinOp.Op = vcat(LinOp.Op,s2)

    push!(LinOp.rhs,0)
  end

  return LinOp.Op,LinOp.rhs
end


"""
  fluidBC!(xy,xyp,dUNodes,dVNodes,dUBCarr,dVBCarr,nNodes,nBCarr,Stiff,F)
    -> Stiff, F

Computes Dirichlet and Neumann boundary conditions
"""
function fluidBC!(mesh::FluidMesh,prob::Problem,
                  LinOp::AbstractLinearOperator)

  fluidDirichlet!(mesh::FluidMesh,prob::Problem,LinOp::LinearOperator)
  LinOp.Op, LinOp.rhs = fluidNeumann(mesh::FluidMesh,prob::Problem,
                                   LinOp::LinearOperator)
end