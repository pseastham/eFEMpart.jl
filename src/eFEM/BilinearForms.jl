# File to be loaded into eFEM

#############################
###### Bilinear Forms #######
#############################

# Mass Matrix
function localMass2D!(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
                       At,phi,dphidx,dphidy,dphids,dphidt,order,parameter...)
  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,nGaussNodes)
  fill!(At,zero(Float64))

  for gpt=1:nGaussNodes
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],order)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,order)
    jac = jacCalc(xN,yN,dphids,dphidt)

    for tj = 1:nNodesPerElm, ti = 1:nNodesPerElm
      integrand = phi[ti]*phi[tj]*jac
      At[ti,tj] += integrand*w[gpt]
    end
  end

  nothing
end

# Laplace2D (no parameter used)
function localLaplace2D!(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
                       At,phi,dphidx,dphidy,dphids,dphidt,order,param)
  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,nGaussNodes)
  fill!(At,zero(Float64))

  for gpt=1:nGaussNodes
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],order)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,order)
    jac = jacCalc(xN,yN,dphids,dphidt)

    for ti = 1:nNodesPerElm, tj = 1:nNodesPerElm
      integrand = (dphidx[ti]*dphidx[tj] + dphidy[ti]*dphidy[tj])/jac
      At[ti,tj] += integrand*w[gpt]
    end
  end

  nothing
end

# Darcy Variable Parameter
function localDarcy2D!(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
  At,phi,dphidx,dphidy,dphids,dphidt,order,param::Vector{T}) where T<:Real
  αNodes = zeros(T,nGaussNodes)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,nGaussNodes)
  fill!(At,zero(T))
  for i=1:nGaussNodes
    αNodes[i] = param[mesh.cm[el].NodeList[i]]
  end

  for gpt=1:nGaussNodes
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],order)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,order)
    jac = jacCalc(xN,yN,dphids,dphidt)
    αg  = shapeEval(αNodes,phi)

    for tj = 1:nNodesPerElm, ti = 1:nNodesPerElm
    integrand = αg*(dphidx[ti]*dphidx[tj] + dphidy[ti]*dphidy[tj])/jac
    At[ti,tj] += integrand*w[gpt]
    end
  end

  nothing
end

# Advection part of Advection-Diffusion Equation
function localAdvDiff2D!(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
                       At,phi,dphidx,dphidy,dphids,dphidt,order,parameter)
  uNodes = zeros(Float64,nGaussNodes)
  vNodes = zeros(Float64,nGaussNodes)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,nGaussNodes)
  fill!(At,zero(Float64))
  for i=1:nNodesPerElm
    uNodes[i] = parameter.u[mesh.cm[el].NodeList[i]]
    vNodes[i] = parameter.v[mesh.cm[el].NodeList[i]]
  end

  for gpt=1:nGaussNodes
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],order)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,order)
    jac = jacCalc(xN,yN,dphids,dphidt)
    ug  = shapeEval(uNodes,phi)
    vg  = shapeEval(vNodes,phi)

    for tj = 1:nNodesPerElm, ti = 1:nNodesPerElm
	    integrand =  ug*phi[ti]*dphidx[tj] + vg*phi[ti]*dphidy[tj]
      At[ti,tj] += integrand*w[gpt]
    end
  end

  nothing
end

# Stokes Constant viscosity
function localStokesConst2D!(mesh,el,xN,yN,w,s,t,Axt,Ayt,Bxt,Byt,
                        phi,dphidx,dphidy,dphids,dphidt,
                        psi,dpsidx,dpsidy,dpsids,dpsidt,
                        parameter...)
  μ0 = parameter[1]                        
  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)

    for ti=1:9
      for tj=1:9
		      integrandX = μ0*dphidx[ti]*dphidx[tj]/jac +
                       μ0*dphidy[ti]*dphidy[tj]/jac
		      Axt[ti,tj] += integrandX*w[gpt]

		      integrandY = μ0*dphidx[ti]*dphidx[tj]/jac +
                       μ0*dphidy[ti]*dphidy[tj]/jac
		      Ayt[ti,tj] += integrandY*w[gpt]
	    end

	    for tk=1:4
	      integrandX = -psi[tk]*dphidx[ti]
	      Bxt[tk,ti] += integrandX*w[gpt]

	      integrandY = -psi[tk]*dphidy[ti]
	      Byt[tk,ti] += integrandY*w[gpt]
	    end
    end
  end
end

# Stokes -- also correct!! -- Variable Viscosity
function localStokesVar2D_WRONG!(mesh,el,xN,yN,w,s,t,Axt,Ayt,A12t,A21t,Bxt,Byt,
                            phi,dphidx,dphidy,dphids,dphidt,
                            psi,dpsidx,dpsidy,dpsids,dpsidt,
                            parameter...)
  μ = parameter[1]

  μNodes = zeros(Float64,9)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    μNodes[i] = μ[mesh.cm[el].NodeList[i]]
  end
  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(A12t,zero(Float64))
  fill!(A21t,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)
    μg  = shapeEval(μNodes,phi)

    for ti=1:9
      for tj=1:9
        integrandX = 2*μg*dphidx[ti]*dphidx[tj]/jac +
                       μg*dphidy[ti]*dphidy[tj]/jac
        Axt[ti,tj] += integrandX*w[gpt]

        integrand12 = μg*dphidx[ti]*dphidy[tj]/jac
        A12t[ti,tj] += integrand12*w[gpt]             

        integrand21 = μg*dphidy[ti]*dphidx[tj]/jac
        A21t[ti,tj] += integrand21*w[gpt]                 

        integrandY =   μg*dphidx[ti]*dphidx[tj]/jac +
                     2*μg*dphidy[ti]*dphidy[tj]/jac
        Ayt[ti,tj] += integrandY*w[gpt]
    end

    for tk=1:4
      integrandX = psi[tk]*dphidx[ti]
      Bxt[tk,ti] -= integrandX*w[gpt]

      integrandY = psi[tk]*dphidy[ti]
      Byt[tk,ti] -= integrandY*w[gpt]
    end
    end
  end

  nothing
end
# Stokes -- correct!! -- Variable Viscosity
function localStokesVar2D!(mesh,el,xN,yN,w,s,t,
                           At,Bt,Ct,Dt,Et,Ft,Gt,Ht,
                           phi,dphidx,dphidy,dphids,dphidt,
                           psi,dpsidx,dpsidy,dpsids,dpsidt,
                           parameter...)
  μ = parameter[1]

  μNodes = zeros(Float64,9)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    μNodes[i] = μ[mesh.cm[el].NodeList[i]]
  end

  fill!(At,zero(Float64))
  fill!(Bt,zero(Float64))
  fill!(Ct,zero(Float64))
  fill!(Et,zero(Float64))
  fill!(Ft,zero(Float64))
  fill!(Gt,zero(Float64))
  fill!(Ht,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)
    μg   = shapeEval(μNodes,phi)

    for ti=1:9
      # vector-valued test function
      for tj=1:9
        # integrand for A
        iA = 2*μg*dphidx[ti]*dphidx[tj]/jac + 
               μg*dphidy[ti]*dphidy[tj]/jac
        At[ti,tj] += iA*w[gpt]

        # integrand for B
        iB = μg*dphidx[ti]*dphidy[tj]/jac
        Bt[ti,tj] += iB*w[gpt]

        # integrand for D -- transpose of B
        iD = μg*dphidx[tj]*dphidy[ti]/jac # iB
        Dt[ti,tj] += iD*w[gpt]        # transpose of B so i and j are switched

        # integrand for E
        iE =   μg*dphidx[ti]*dphidx[tj]/jac + 
             2*μg*dphidy[ti]*dphidy[tj]/jac
        Et[ti,tj] += iE*w[gpt] 
      end

      # scalar-valued test function
      for tj=1:4
        # integrand for C
        iC = psi[tj]*dphidx[ti] # C is Nu x Np block matrix, so i and j are correct
        Ct[ti,tj] -= iC*w[gpt]

        # integrand for F
        iF = psi[tj]*dphidy[ti] # F is Nu x Np block matrix, so i and j are correct
        Ft[ti,tj] -= iF*w[gpt]

        # integrand for G -- transpose of C
        iG = iC
        Gt[tj,ti] -= iG*w[gpt]    # G is Np x Nu block matrix, so i and j are switched

        # integrand for H -- transpose of H
        iH = iF
        Ht[tj,ti] -= iH*w[gpt]    # H is Np x Nu block matrix, so i and j are switched
      end
    end
  end

  nothing
end
# Stokes -- correct!! -- Variable Viscosity
function localStokesAlternateVar2D!(mesh,el,xN,yN,w,s,t,
                                    At,Bt,Ct,Dt,Et,Ft,Gt,Ht,
                                    phi,dphidx,dphidy,dphids,dphidt,
                                    psi,dpsidx,dpsidy,dpsids,dpsidt,
                                    parameter...)
  μ = parameter[1]

  μNodes = zeros(Float64,9)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    μNodes[i] = μ[mesh.cm[el].NodeList[i]]
  end

  fill!(At,zero(Float64))
  fill!(Bt,zero(Float64))
  fill!(Ct,zero(Float64))
  fill!(Et,zero(Float64))
  fill!(Ft,zero(Float64))
  fill!(Gt,zero(Float64))
  fill!(Ht,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)
    
    μg    = shapeEval(μNodes,phi)
    dμgdx = shapeEval(μNodes,dphidx)
    dμgdy = shapeEval(μNodes,dphidy)

    for ti=1:9
      # vector-valued test function
      for tj=1:9
        # integrand for A
        iA = μg*(dphidx[ti]*dphidx[tj] + dphidy[ti]*dphidy[tj])/jac -
             phi[ti]*(dμgdx*dphidx[tj] + dμgdy*dphidy[tj])/jac
        At[ti,tj] += iA*w[gpt]

        # integrand for B
        iB = 0.0
        Bt[ti,tj] += iB*w[gpt]

        # integrand for D -- transpose of B
        iD = 0.0
        Dt[ti,tj] += iD*w[gpt]

        # integrand for E
        iE = iA
        Et[ti,tj] += iE*w[gpt] 
      end

      # scalar-valued test function
      for tj=1:4
        # integrand for C
        iC = -psi[tj]*dphidx[ti] # C is Nu x Np block matrix, so i and j are correct
        Ct[ti,tj] += iC*w[gpt]

        # integrand for F
        iF = -psi[tj]*dphidy[ti] # F is Nu x Np block matrix, so i and j are correct
        Ft[ti,tj] += iF*w[gpt]

        # integrand for G -- transpose of C
        iG = iC
        Gt[tj,ti] += iG*w[gpt]    # G is Np x Nu block matrix, so i and j are switched

        # integrand for H -- transpose of H
        iH = iF
        Ht[tj,ti] += iH*w[gpt]    # H is Np x Nu block matrix, so i and j are switched
      end
    end
  end

nothing
end

# Brinkman Constant parameter
function localBrinkmanConst2D!(mesh,el,xN,yN,w,s,t,Axt,Ayt,Bxt,Byt,
                               phi,dphidx,dphidy,dphids,dphidt,
                               psi,dpsidx,dpsidy,dpsids,dpsidt,
                               parameter...)
  μ = parameter[1]
  α = parameter[2]

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)

    for ti=1:9
      for tj=1:9
        integrandX = μ*dphidx[ti]*dphidx[tj]/jac +
                     μ*dphidy[ti]*dphidy[tj]/jac +
                     α*phi[ti]*phi[tj]*jac
        Axt[ti,tj] += integrandX*w[gpt]

        integrandY = μ*dphidx[ti]*dphidx[tj]/jac +
                     μ*dphidy[ti]*dphidy[tj]/jac +
                     α*phi[ti]*phi[tj]*jac
        Ayt[ti,tj] += integrandY*w[gpt]
      end

      for tk=1:4
        integrandX = psi[tk]*dphidx[ti]
        Bxt[tk,ti] -= integrandX*w[gpt]

        integrandY = psi[tk]*dphidy[ti]
        Byt[tk,ti] -= integrandY*w[gpt]
      end
    end
  end
end


# Brinkman Variable parameter
function localBrinkmanVar2D!(mesh,el,xN,yN,w,s,t,Axt,Ayt,A12t,A21t,Bxt,Byt,
                             phi,dphidx,dphidy,dphids,dphidt,
                             psi,dpsidx,dpsidy,dpsids,dpsidt,
                             parameter...)
  μ = parameter[1]
  α = parameter[2]

  μNodes = zeros(Float64,9)
  αNodes  = zeros(Float64,9)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    μNodes[i] = μ[mesh.cm[el].NodeList[i]]
    αNodes[i] = α[mesh.cm[el].NodeList[i]]
  end
  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(A12t,zero(Float64))
  fill!(A21t,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)

    μg = shapeEval(μNodes,phi)
    αg = shapeEval(αNodes,phi)

    for ti=1:9
      for tj=1:9
        integrandX  = 2*μg*dphidx[ti]*dphidx[tj]/jac +
                        μg*dphidy[ti]*dphidy[tj]/jac +
                        αg*phi[ti]*phi[tj]*jac
        Axt[ti,tj] += integrandX*w[gpt]

        integrand12  = μg*dphidx[tj]*dphidy[ti]/jac
        A12t[ti,tj] += integrand12*w[gpt]

        integrand21  = μg*dphidy[tj]*dphidx[ti]/jac
        A21t[ti,tj] += integrand21*w[gpt]

        integrandY   =   μg*dphidx[ti]*dphidx[tj]/jac +
                       2*μg*dphidy[ti]*dphidy[tj]/jac +
                         αg*phi[ti]*phi[tj]*jac
        Ayt[ti,tj]  += integrandY*w[gpt]
      end

      for tk=1:4
        integrandX = psi[tk]*dphidx[ti]
        Bxt[tk,ti] -= integrandX*w[gpt]

        integrandY = psi[tk]*dphidy[ti]
        Byt[tk,ti] -= integrandY*w[gpt]
      end
    end
  end
end

# Brinkman Multiphase 2D Constant Parameter
function localBrinkmanMPConst2D!(mesh,el,xN,yN,w,s,t,Axt,Ayt,Bxt,Byt,
                                  phi,dphidx,dphidy,dphids,dphidt,
                                  psi,dpsidx,dpsidy,dpsids,dpsidt,
                                  parameter...)
  α1 = parameter[1]
  α2 = parameter[2]
  α3 = parameter[3]

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)

  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)

    for ti=1:9
      for tj=1:9
        integrandX = α1*dphidx[ti]*dphidx[tj]/jac +
                     α1*dphidy[ti]*dphidy[tj]/jac +
                     α2*phi[ti]*phi[tj]*jac
        Axt[ti,tj] += integrandX*w[gpt]

        integrandY = α1*dphidx[ti]*dphidx[tj]/jac +
                     α1*dphidy[ti]*dphidy[tj]/jac +
                     α2*phi[ti]*phi[tj]*jac
        Ayt[ti,tj] += integrandY*w[gpt]
      end

      for tk=1:4
        integrandX = α3*psi[tk]*dphidx[ti]
        Bxt[tk,ti] -= integrandX*w[gpt]

        integrandY = α3*psi[tk]*dphidy[ti]
        Byt[tk,ti] -= integrandY*w[gpt]
      end
    end
  end
end

# Brinkman Multiphase -- wrong!! -- 2D Variable Parameter
function localBrinkmanMPVar2D!(mesh,el,xN,yN,w,s,t,Axt,Ayt,A12t,A21t,Bxt,Byt,
                               phi,dphidx,dphidy,dphids,dphidt,
                               psi,dpsidx,dpsidy,dpsids,dpsidt,
                               parameter...)
  α1 = parameter[1]
  α2 = parameter[2]
  α3 = parameter[3]

  α1Nodes = zeros(Float64,9)
  α2Nodes = zeros(Float64,9)
  α3Nodes = zeros(Float64,9)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    α1Nodes[i] = α1[mesh.cm[el].NodeList[i]]
    α2Nodes[i] = α2[mesh.cm[el].NodeList[i]]
    α3Nodes[i] = α3[mesh.cm[el].NodeList[i]]
  end

  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(A12t,zero(Float64))
  fill!(A21t,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

    for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)

    α1g = shapeEval(α1Nodes,phi)
    α2g = shapeEval(α2Nodes,phi)
    α3g = shapeEval(α3Nodes,phi)
    yg  = shapeEval(yN,phi)

    for ti=1:9
      for tj=1:9
        integrandX = 2*α1g*dphidx[ti]*dphidx[tj]/jac +
                       α1g*dphidy[ti]*dphidy[tj]/jac +
                       α2g*phi[ti]*phi[tj]*jac
        Axt[ti,tj] += integrandX*w[gpt]

        integrand12  = α1g*dphidx[tj]*dphidy[ti]/jac
        A12t[ti,tj] += integrand12*w[gpt]
        
        integrand21  = α1g*dphidy[tj]*dphidx[ti]/jac
        A21t[ti,tj] += integrand21*w[gpt]
        
        integrandY =   α1g*dphidx[ti]*dphidx[tj]/jac +
                     2*α1g*dphidy[ti]*dphidy[tj]/jac +
                       α2g*phi[ti]*phi[tj]*jac
        Ayt[ti,tj] += integrandY*w[gpt]
      end

      for tk=1:4
        integrandX = α3g*psi[tk]*dphidx[ti]
        Bxt[tk,ti] -= integrandX*w[gpt]

        integrandY = α3g*psi[tk]*dphidy[ti]
        Byt[tk,ti] -= integrandY*w[gpt]
      end
    end
  end
end
# Multiphase Brinkman (MPB) -- correct!! -- 2D Variable Parameter
function localMPBVar2D!(mesh,el,xN,yN,w,s,t,
                        At,Ct,Et,Ft,Gt,Ht,
                        phi,dphidx,dphidy,dphids,dphidt,
                        psi,dpsidx,dpsidy,dpsids,dpsidt,
                        parameter...)
  # make parameters easier to deal with
  α1 = parameter[1]
  α2 = parameter[2]
  α3 = parameter[3]

  # temporary arrays for variable-coefficients
  α1Nodes = zeros(Float64,9)
  α2Nodes = zeros(Float64,9)
  α3Nodes = zeros(Float64,9)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    α1Nodes[i] = α1[mesh.cm[el].NodeList[i]]
    α2Nodes[i] = α2[mesh.cm[el].NodeList[i]]
    α3Nodes[i] = α3[mesh.cm[el].NodeList[i]]
  end

  fill!(At,zero(Float64))
  fill!(Ct,zero(Float64))
  fill!(Et,zero(Float64))
  fill!(Ft,zero(Float64))
  fill!(Gt,zero(Float64))
  fill!(Ht,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)

    α1g = shapeEval(α1Nodes,phi)
    α2g = shapeEval(α2Nodes,phi)
    α3g = shapeEval(α3Nodes,phi)
    
    dα1gdx = shapeEval(α1Nodes,dphidx)
    dα1gdy = shapeEval(α1Nodes,dphidy)
    dα3gdx = shapeEval(α3Nodes,dphidx)
    dα3gdy = shapeEval(α3Nodes,dphidy)

    for ti=1:9
      # vector-valued test function
      for tj=1:9
        # integrand for A
        iA = α1g*(dphidx[ti]*dphidx[tj] + dphidy[ti]*dphidy[tj])/jac + 
             phi[ti]*(dα1gdx*dphidx[tj] + dα1gdy*dphidy[tj])/jac +
             α2g*phi[ti]*phi[tj]*jac
        At[ti,tj] += iA*w[gpt]

        # integrand for E -- same as A
        iE = iA
        Et[ti,tj] += iE*w[gpt]    # A and E are same block matrices
      end

      # scalar-valued test function
      for tj=1:4
        # integrand for C
        iC = -psi[tj]*(α3g*dphidx[ti] + dα3gdx*phi[ti]) # C is Nu x Np block matrix, so i and j are correct
        Ct[ti,tj] += iC*w[gpt]

        # integrand for F
        iF = -psi[tj]*(α3g*dphidy[ti] + dα3gdy*phi[ti]) # F is Nu x Np block matrix, so i and j are correct
        Ft[ti,tj] += iF*w[gpt]

        # integrand for G
        iG = -psi[tj]*dphidx[ti]
        Gt[tj,ti] += iG*w[gpt]    # G is Np x Nu block matrix, so i and j are switched

        # integrand for H
        iH = -psi[tj]*dphidy[ti]
        Ht[tj,ti] += iH*w[gpt]    # H is Np x Nu block matrix, so i and j are switched
      end
    end
  end
end

# Constant diffusion Laplace Axisymmetric
function localLaplaceConstAS!(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
                         At,phi,dphidx,dphidy,dphids,dphidt,order,parameter...)
  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,nGaussNodes)
  # zero out At
  fill!(At,zero(Float64))

  for gpt=1:nGaussNodes
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],order)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,order)
    jac = jacCalc(xN,yN,dphids,dphidt)
    yg = shapeEval(yN,phi)

    for ti = 1:nNodesPerElm, tj = 1:nNodesPerElm
      integrand = (dphidx[ti]*dphidx[tj] + dphidy[ti]*dphidy[tj])/jac
      At[ti,tj] += integrand*w[gpt]*yg
    end
  end

  nothing
end

# Constant diffusion Laplace Axisymmetric variable coefficient
function localLaplaceVarAS!(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
                            At,phi,dphidx,dphidy,dphids,dphidt,order,parameter)
  
  κNodes = zeros(Float64,nGaussNodes)                          
  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,nGaussNodes)
  # zero out At
  fill!(At,zero(Float64))
  for i=1:nGaussNodes
    κNodes[i] = parameter[mesh.cm[el].NodeList[i]]
  end
  for gpt=1:nGaussNodes
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],order)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,order)
    jac = jacCalc(xN,yN,dphids,dphidt)
    κg = shapeEval(κNodes,phi)
    yg = shapeEval(yN,phi)

    for ti = 1:nNodesPerElm, tj = 1:nNodesPerElm
      integrand = κg*(dphidx[ti]*dphidx[tj] + dphidy[ti]*dphidy[tj])/jac
      At[ti,tj] += integrand*w[gpt]*yg
    end
  end

  nothing
end

function localAdvDiffAS!(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
                         At,phi,dphidx,dphidy,dphids,dphidt,order,parameter)
  U = parameter.u
  V = parameter.v

  UNodes = zeros(Float64,nGaussNodes)
  VNodes = zeros(Float64,nGaussNodes)
  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,nGaussNodes)
  fill!(At,zero(Float64))

  for i=1:nGaussNodes
    UNodes[i] = U[mesh.cm[el].NodeList[i]]
    VNodes[i] = V[mesh.cm[el].NodeList[i]]
  end

  for gpt=1:nGaussNodes
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],order)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,order)
    jac = jacCalc(xN,yN,dphids,dphidt)
    ug  = shapeEval(UNodes,phi)
    vg  = shapeEval(VNodes,phi)
    yg  = shapeEval(yN,phi)

    for ti = 1:nNodesPerElm, tj = 1:nNodesPerElm
      integrand =  ug*phi[ti]*dphidx[tj] + vg*phi[ti]*dphidy[tj]
      At[ti,tj] += integrand*w[gpt]*yg
    end
  end

  nothing
end

function localStokesConstAS!(mesh,el,xN,yN,w,s,t,Axt,Ayt,Bxt,Byt,
                             phi,dphidx,dphidy,dphids,dphidt,
                             psi,dpsidx,dpsidy,dpsids,dpsidt,
                             parameter...)
  μ0 = parameter[1]
  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)

  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:9
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)
    yg  = shapeEval(yN,phi)

    for ti=1:9
      for tj=1:9
        integrandX = μ0*dphidx[ti]*dphidx[tj]*yg/jac +
                     μ0*dphidy[ti]*dphidy[tj]*yg/jac
        Axt[ti,tj] += integrandX*w[gpt]

        integrandY = μ0*dphidx[ti]*dphidx[tj]*yg/jac +
                     μ0*dphidy[ti]*dphidy[tj]*yg/jac +
                     μ0*phi[ti]*phi[tj]*jac/yg
        Ayt[ti,tj] += integrandY*w[gpt]
      end

      for tk=1:4
        integrandX = psi[tk]*dphidx[ti]*yg
        Bxt[tk,ti] -= integrandX*w[gpt]

        integrandY = psi[tk]*dphidy[ti]*yg + psi[tk]*phi[ti]*jac
        Byt[tk,ti] -= integrandY*w[gpt]
      end
    end
  end
end

function localStokesVarAS!(mesh,el,xN,yN,w,s,t,Axt,Ayt,A12t,A21t,Bxt,Byt,
                            phi,dphidx,dphidy,dphids,dphidt,
                            psi,dpsidx,dpsidy,dpsids,dpsidt,
                            parameter...)
  μ = parameter[1]
  μNodes = zeros(Float64,9)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    μNodes[i] = μ[mesh.cm[el].NodeList[i]]
  end

  nGpt = length(s)

  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(A12t,zero(Float64))
  fill!(A21t,zero(Float64))  
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:nGpt
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)
    yg  = shapeEval(yN,phi)
    μg  = shapeEval(μNodes,phi)

    for ti=1:9
      for tj=1:9
        integrandX = 2.0*μg*dphidx[ti]*dphidx[tj]*yg/jac +
                         μg*dphidy[ti]*dphidy[tj]*yg/jac
        Axt[ti,tj] += integrandX*w[gpt]

        # NOTE: I'm pretty sure the A21t is wrong, but there's something later on in the assembly process 
        #       that fixes this bug and so for the love of god I know it's wrong but please don't change it
        #       because it's fixed later on by another bug...will have to go through and check carefully 
        #       at a later date.
        # 
        # integrand21 should be = μg*dphidx[tj]*dphidy[ti]*yg/jac. Note indices
        integrand12 = μg*dphidx[tj]*dphidy[ti]*yg/jac   
        A12t[ti,tj] += integrand12*w[gpt]

        integrand21 = μg*dphidx[ti]*dphidy[tj]*yg/jac
        A21t[ti,tj] += integrand21*w[gpt]

        integrandY =     μg*dphidx[ti]*dphidx[tj]*yg/jac +
                     2.0*μg*dphidy[ti]*dphidy[tj]*yg/jac +
                     2.0*μg*phi[ti]*phi[tj]*jac/yg
        Ayt[ti,tj] += integrandY*w[gpt]
      end

      for tk=1:4
        integrandX = psi[tk]*dphidx[ti]*yg
        Bxt[tk,ti] -= integrandX*w[gpt]

        integrandY = psi[tk]*dphidy[ti]*yg + psi[tk]*phi[ti]*jac
        Byt[tk,ti] -= integrandY*w[gpt]
      end
    end
  end
end
function localStokesVarASAlternate!(mesh,el,xN,yN,w,s,t,Axt,Ayt,Bxt,Byt,
                                    phi,dphidx,dphidy,dphids,dphidt,
                                    psi,dpsidx,dpsidy,dpsids,dpsidt,
                                    parameter...)
  μ = parameter[1]
  μNodes = zeros(Float64,9)

  nGpt = length(s)

  # generate local stiffness matrices
  getNodes!(xN,yN,mesh,el,9)
  for i=1:9
    μNodes[i] = μ[mesh.cm[el].NodeList[i]]
  end

  fill!(Axt,zero(Float64))
  fill!(Ayt,zero(Float64))
  fill!(Bxt,zero(Float64))
  fill!(Byt,zero(Float64))

  for gpt=1:nGpt
    shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],2)
    shape2D!(psi,dpsids,dpsidt,s[gpt],t[gpt],1)
    derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xN,yN,2)
    derivShape2D!(psi,dpsidx,dpsidy,dpsids,dpsidt,s[gpt],t[gpt],xN[1:4],yN[1:4],1)
    jac = jacCalc(xN,yN,dphids,dphidt)
    yg  = shapeEval(yN,phi)
    μg  = shapeEval(μNodes,phi)

    dμgdx = shapeEval(μNodes,dphidx)
    dμgdy = shapeEval(μNodes,dphidy)

    for ti=1:9
      for tj=1:9
        integrandX = μg*dphidx[ti]*dphidx[tj]*yg/jac +
                     μg*dphidy[ti]*dphidy[tj]*yg/jac -
                     dμgdx*phi[tj]*dphidx[ti]*yg/jac -
                     dμgdy*phi[tj]*dphidy[ti]*yg/jac
        Axt[ti,tj] += integrandX*w[gpt]

        integrandY = integrandX + μg*phi[ti]*phi[tj]/yg*jac 
        Ayt[ti,tj] += integrandY*w[gpt]
      end

      for tk=1:4
        integrandX = psi[tk]*dphidx[ti]*yg
        Bxt[tk,ti] -= integrandX*w[gpt]

        integrandY = psi[tk]*dphidy[ti]*yg + psi[tk]*phi[ti]*jac
        Byt[tk,ti] -= integrandY*w[gpt]
      end
    end
  end
end

################################
###### Support Functions #######
################################

function getNodes(mesh,el,nGaussNodes)
  xNodes = Array{Float64}(nGaussNodes); yNodes = Array{Float64}(nGaussNodes)
  for i=1:nGaussNodes
    xNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].x
    yNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].y
  end
  return xNodes,yNodes
end

function getNodes!(xNodes,yNodes,mesh,el,nGaussNodes)
  for i=1:nGaussNodes
    xNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].x
    yNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].y
  end
end