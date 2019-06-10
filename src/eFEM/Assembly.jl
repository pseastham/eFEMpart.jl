# to be loaded into eFEMpart
# Assembles global linear system from local element matrices
# which are represented in functional form by 'localmat' argument

########################
###### Assembler #######
########################

"""
  assembleScalar(mesh,localmat,parameter)

Assembles sparse matrix for scalar PDEs

ASSUMPTIONS:
- mesh is either Mesh or FluidMesh type
- localmat is functional form of PDE for local element
- parameter can be either number or array (variable coefficient problem)

INPUT: 
  mesh:       [AbstractMesh]
                  mesh type that PDE is being solved on
  localmat:   [Function]
                  function dictating what PDE is being solved 
  parameter:  [Real or Vector{Real}]
                  parameter for PDE, can be either of type Real (constant coefficient 
                    PDE) or Vector{Real} (variable coefficient PDE)

OUTPUT:
  -           [SparseArray]
                  sparse array for linear system to be solved
"""
function assembleScalar(mesh,localmat,parameter)
  # derive parameters
  nNodes = length(mesh.xy)
  nElm   = length(mesh.cm)
  (mesh.order==:Linear ? order=1 : order=2)
  (order==1 ? nNodesPerElm=4 : nNodesPerElm=9)
  nGaussNodes = nNodesPerElm  # possibly not true for axisymmetric quadrature

  #tempArr = Array{Float64}(nNodesPerElm)
  xN      = Array{Float64}(undef,nNodesPerElm)
  yN      = Array{Float64}(undef,nNodesPerElm)
  phi     = Array{Float64}(undef,nNodesPerElm)
  dphidx  = Array{Float64}(undef,nNodesPerElm)
  dphidy  = Array{Float64}(undef,nNodesPerElm)
  dphids  = Array{Float64}(undef,nNodesPerElm)
  dphidt  = Array{Float64}(undef,nNodesPerElm)

  At  = Array{Float64}(undef,nNodesPerElm,nNodesPerElm)

  I = Array{Int64}(undef,nElm,nNodesPerElm,nNodesPerElm)
  J = Array{Int64}(undef,nElm,nNodesPerElm,nNodesPerElm)
  S = Array{Float64}(undef,nElm,nNodesPerElm,nNodesPerElm)

  # compute phi,dphids,dphidt
  w,s,t = GaussQuadPoints2D(order+1)

  # construct local matrices
  for el = 1:nElm
    localmat(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
             At,phi,dphidx,dphidy,dphids,dphidt,order,parameter)
    for ti = 1:nNodesPerElm, tj = 1:nNodesPerElm
      I[el,ti,tj] = mesh.cm[el].NodeList[ti]
      J[el,ti,tj] = mesh.cm[el].NodeList[tj]
      S[el,ti,tj] = At[ti,tj]
    end
  end

  # assemble global matrix
  return sparse(vec(I),vec(J),vec(S),nNodes,nNodes)
end
# implements dirichlet boundary condition information before neumann
# and affects F simultaneously
# as per https://www.math.colostate.edu/~bangerth/videos.676.21.6.html
function assembleScalarandRHS_2b(mesh,localmat,parameter,prob)
  # derive parameters
  nNodes = length(mesh.xy)
  nElm   = length(mesh.cm)
  (mesh.order==:Linear ? order=1 : order=2)
  (order==1 ? nNodesPerElm=4 : nNodesPerElm=9)
  nGaussNodes = nNodesPerElm  # possibly not true for axisymmetric quadrature

  #tempArr = Array{Float64}(nNodesPerElm)
  xN      = Array{Float64}(undef,nNodesPerElm)
  yN      = Array{Float64}(undef,nNodesPerElm)
  phi     = Array{Float64}(undef,nNodesPerElm)
  dphidx  = Array{Float64}(undef,nNodesPerElm)
  dphidy  = Array{Float64}(undef,nNodesPerElm)
  dphids  = Array{Float64}(undef,nNodesPerElm)
  dphidt  = Array{Float64}(undef,nNodesPerElm)

  At  = Array{Float64}(undef,nNodesPerElm,nNodesPerElm)

  # indice arrays preallocation
  I = Array{Int64}(undef,nElm,nNodesPerElm,nNodesPerElm)
  J = Array{Int64}(undef,nElm,nNodesPerElm,nNodesPerElm)
  S = Array{Float64}(undef,nElm,nNodesPerElm,nNodesPerElm)

  # forcing array preallocation
  F = zeros(Float64,nNodes)

  # for indexing how many times boundary is accessed
  bcCOUNT = zeros(Float64,length(prob.bcid[:dNodes]))

  # compute phi,dphids,dphidt
  w,s,t = GaussQuadPoints2D(order+1)

  # construct local matrices
  # HERE IS WHERE YOU IMPLEMENT BOUNDARY DATA
  for el = 1:nElm
    localmat(mesh,el,xN,yN,w,s,t,nGaussNodes,nNodesPerElm,
             At,phi,dphidx,dphidy,dphids,dphidt,order,parameter)
    # local assembly -- nothing changes
    for ti = 1:nNodesPerElm, tj = 1:nNodesPerElm
      I[el,ti,tj] = mesh.cm[el].NodeList[ti]
      J[el,ti,tj] = mesh.cm[el].NodeList[tj]
      S[el,ti,tj] = At[ti,tj]
    end

    # local F assembly
    c = mesh.cm[el].NodeList
    for gpt=1:nNodesPerElm
      phi,_,_,jac = derivShape2D(s[gpt],t[gpt],xN,yN,order)
	    for i=1:nNodesPerElm
	      fcoefs      = prob.bcval[:forcing][mesh.cm[el].NodeList]
	      fgpt        = shapeEval(fcoefs,phi)
	      F[mesh.cm[el].NodeList[i]] += fgpt*phi[i]*w[gpt]*jac
	    end
	  end

    # check for boundary data & modify local stiffness and forcing arrays accordingly
    for ti=1:nNodesPerElm
      if mesh.cm[el].NodeList[ti] in prob.bcid[:dNodes]
        #find index of bcid
        bcIND = findfirst(x->x==mesh.cm[el].NodeList[ti], prob.bcid[:dNodes])
        bcCOUNT[bcIND] += 1

        # subtract matvec from F
        for tj=1:nNodesPerElm
          F[mesh.cm[el].NodeList[tj]] -= S[el,tj,ti]*prob.bcval[:dBC][bcIND]
          S[el,tj,ti] = 0.0
          S[el,ti,tj] = 0.0
        end

        # place 1 in diagonal of local stiffness
        S[el,ti,ti] = 1.0

        # replace Frocing with dirichlet value
        F[mesh.cm[el].NodeList[ti]] = prob.bcval[:dBC][bcIND]
      end
    end
  end

  F[prob.bcid[:dNodes]] .*= bcCOUNT

  # assemble global matrix
  return sparse(vec(I),vec(J),vec(S),nNodes,nNodes), F
end

"""
  assembleHalfFluid(mesh,localmat,parameter...)

Assembles sparse matrix for Fluids with no cross-terms

ASSUMPTIONS:
- mesh is either Mesh or FluidMesh type
- localmat is functional form of PDE for local element
- parameter can be either number or array (variable coefficient problem),
  also can be a tuple of parameters depending on expectation of localmat

INPUT: 
  mesh:       [AbstractMesh]
                  mesh type that PDE is being solved on
  localmat:   [Function]
                  function dictating what PDE is being solved 
  parameter:  [Real or Vector{Real}]
                  parameter for PDE, can be either of type Real (constant coefficient 
                    PDE) or Vector{Real} (variable coefficient PDE)

OUTPUT:
  -           [SparseArray]
                  sparse array for linear system to be solved
"""
function assembleHalfFluid(mesh,localmat,parameter...)
  nElm = length(mesh.cm)
  nUNodes = length(mesh.xy); nPNodes = length(mesh.xyp)
  nGaussNodes = 9
  nNodesPerElmQ1 = 4
  nNodesPerElmQ2 = 9

  xN = Array{Float64}(undef,nNodesPerElmQ2)
  yN = Array{Float64}(undef,nNodesPerElmQ2)
  phi    = Array{Float64}(undef,nNodesPerElmQ2)
  dphids = Array{Float64}(undef,nNodesPerElmQ2)
  dphidt = Array{Float64}(undef,nNodesPerElmQ2) 
  dphidx = Array{Float64}(undef,nNodesPerElmQ2)
  dphidy = Array{Float64}(undef,nNodesPerElmQ2)
  psi    = Array{Float64}(undef,nNodesPerElmQ1)
  dpsids = Array{Float64}(undef,nNodesPerElmQ1)
  dpsidt = Array{Float64}(undef,nNodesPerElmQ1) 
  dpsidx = Array{Float64}(undef,nNodesPerElmQ1)
  dpsidy = Array{Float64}(undef,nNodesPerElmQ1)

  Axt = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Ayt = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Bxt = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)
  Byt = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)

  IA1 = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  JA1 = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  SA1 = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  IA2 = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  JA2 = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  SA2 = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  IBx1 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JBx1 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  SBx1 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  IBy1 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JBy1 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  SBy1 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  IBx2 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JBx2 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  SBx2 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  IBy2 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JBy2 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  SBy2 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)

  # for axisymmetric, better to use 4x4 gauss?
  w,s,t = GaussQuadPoints2D(3)

  for el = 1:nElm
    localmat(mesh,el,xN,yN,w,s,t,Axt,Ayt,Bxt,Byt,
              phi,dphidx,dphidy,dphids,dphidt,
              psi,dpsidx,dpsidy,dpsids,dpsidt,
              parameter...)

    for ti = 1:nNodesPerElmQ2, tj = 1:nNodesPerElmQ2
      # Axt
      IA1[el,ti,tj] = mesh.cm[el].NodeList[ti]
      JA1[el,ti,tj] = mesh.cm[el].NodeList[tj]
      SA1[el,ti,tj] = Axt[ti,tj]
      # Ayt
      IA2[el,ti,tj] = mesh.cm[el].NodeList[ti] + nUNodes
      JA2[el,ti,tj] = mesh.cm[el].NodeList[tj] + nUNodes
      SA2[el,ti,tj] = Ayt[ti,tj]
    end
    for tk = 1:nNodesPerElmQ1, ti=1:nNodesPerElmQ2
      # Bx
      IBx1[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      JBx1[el,tk,ti] = mesh.cm[el].NodeList[ti]
      SBx1[el,tk,ti] = Bxt[tk,ti]
      # By
      IBy1[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      JBy1[el,tk,ti] = mesh.cm[el].NodeList[ti]  + nUNodes
      SBy1[el,tk,ti] = Byt[tk,ti]
      # Bx'
      IBx2[el,tk,ti] = mesh.cm[el].NodeList[ti]
      JBx2[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      SBx2[el,tk,ti] = Bxt[tk,ti]
      # By'
      IBy2[el,tk,ti] = mesh.cm[el].NodeList[ti]  + nUNodes
      JBy2[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      SBy2[el,tk,ti] = Byt[tk,ti]
    end
  end

  I = vcat(vec(IA1),vec(IA2),vec(IBx1),vec(IBx2),vec(IBy1),vec(IBy2))
  J = vcat(vec(JA1),vec(JA2),vec(JBx1),vec(JBx2),vec(JBy1),vec(JBy2))
  S = vcat(vec(SA1),vec(SA2),vec(SBx1),vec(SBx2),vec(SBy1),vec(SBy2))

  return sparse(I,J,S,2*nUNodes+nPNodes,2*nUNodes+nPNodes)
end
"""
  assembleFullFluid_WRONG(mesh,localmat,parameter...)

actually not wrong... ?????
"""
function assembleFullFluid_WRONG(mesh,localmat,parameter...)
  nElm = length(mesh.cm)
  nUNodes = length(mesh.xy); nPNodes = length(mesh.xyp)
  nGaussNodes = 9
  nNodesPerElmQ1 = 4
  nNodesPerElmQ2 = 9

  w,s,t = GaussQuadPoints2D(3)

  tempArr4 = Array{Float64}(undef,nNodesPerElmQ1)
  tempArr9 = Array{Float64}(undef,nNodesPerElmQ2)

  xN     = copy(tempArr9); yN     = copy(tempArr9)

  phi    = copy(tempArr9); dphids = copy(tempArr9); dphidt=copy(tempArr9)
  dphidx = copy(tempArr9); dphidy = copy(tempArr9)

  psi    = copy(tempArr4); dpsids = copy(tempArr4); dpsidt = copy(tempArr4)
  dpsidx = copy(tempArr4); dpsidy = copy(tempArr4)

  # top left
  Axt  = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  # top right
  A12t = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  # bottom left
  A21t = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  # bottom right
  Ayt  = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Bxt = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)
  Byt = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)

  IA11 = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  JA11 = copy(IA11); SA11 = copy(IA11);
  IA12 = copy(IA11); JA12 = copy(IA11); SA12 = copy(IA11)
  IA21 = copy(IA11); JA21 = copy(IA11); SA21 = copy(IA11)
  IA22 = copy(IA11); JA22 = copy(IA11); SA22 = copy(IA11)

  IBx1 = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JBx1 = copy(IBx1); SBx1 = copy(IBx1)
  IBy1 = copy(IBx1); JBy1 = copy(IBx1); SBy1 = copy(IBx1)
  IBx2 = copy(IBx1); JBx2 = copy(IBx1); SBx2 = copy(IBx1)
  IBy2 = copy(IBx1); JBy2 = copy(IBx1); SBy2 = copy(IBx1)

  for el = 1:nElm
    localmat(mesh,el,xN,yN,w,s,t,Axt,Ayt,A12t,A21t,Bxt,Byt,
              phi,dphidx,dphidy,dphids,dphidt,
              psi,dpsidx,dpsidy,dpsids,dpsidt,
              parameter...)

    for ti = 1:nNodesPerElmQ2, tj = 1:nNodesPerElmQ2
      # Axt
      IA11[el,ti,tj] = mesh.cm[el].NodeList[ti]
      JA11[el,ti,tj] = mesh.cm[el].NodeList[tj]
      SA11[el,ti,tj] = Axt[ti,tj]
      # A12t
      IA12[el,ti,tj] = mesh.cm[el].NodeList[ti]
      JA12[el,ti,tj] = mesh.cm[el].NodeList[tj] + nUNodes
      SA12[el,ti,tj] = A12t[ti,tj]
      # A21t
      IA21[el,ti,tj] = mesh.cm[el].NodeList[ti] + nUNodes
      JA21[el,ti,tj] = mesh.cm[el].NodeList[tj]
      SA21[el,ti,tj] = A21t[ti,tj]
      # Ayt
      IA22[el,ti,tj] = mesh.cm[el].NodeList[ti] + nUNodes
      JA22[el,ti,tj] = mesh.cm[el].NodeList[tj] + nUNodes
      SA22[el,ti,tj] = Ayt[ti,tj]
    end
    for tk = 1:nNodesPerElmQ1, ti=1:nNodesPerElmQ2
      # Bx
      IBx1[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      JBx1[el,tk,ti] = mesh.cm[el].NodeList[ti]
      SBx1[el,tk,ti] = Bxt[tk,ti]
      # By
      IBy1[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      JBy1[el,tk,ti] = mesh.cm[el].NodeList[ti]  + nUNodes
      SBy1[el,tk,ti] = Byt[tk,ti]
      # Bx'
      IBx2[el,tk,ti] = mesh.cm[el].NodeList[ti]
      JBx2[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      SBx2[el,tk,ti] = Bxt[tk,ti]
      # By'
      IBy2[el,tk,ti] = mesh.cm[el].NodeList[ti]  + nUNodes
      JBy2[el,tk,ti] = mesh.cmp[el].NodeList[tk] + 2*nUNodes
      SBy2[el,tk,ti] = Byt[tk,ti]
    end
  end

  I = vcat(vec(IA11),vec(IA12),vec(IA21),vec(IA22),vec(IBx1),vec(IBx2),vec(IBy1),vec(IBy2))
  J = vcat(vec(JA11),vec(JA12),vec(JA21),vec(JA22),vec(JBx1),vec(JBx2),vec(JBy1),vec(JBy2))
  S = vcat(vec(SA11),vec(SA12),vec(SA21),vec(SA22),vec(SBx1),vec(SBx2),vec(SBy1),vec(SBy2))

  return sparse(I,J,S,2*nUNodes+nPNodes,2*nUNodes+nPNodes)
end
"""
  assembleFullFluid(mesh,localmat,parameter...)

Assembles full 3x3 block matrix for fluid problems
"""
# Assembles full 3x3 block matrix for fluid problems
function assembleFullFluid(mesh,localmat,parameter...)
  nElm = length(mesh.cm)
  nUNodes = length(mesh.xy); nPNodes = length(mesh.xyp)
  nGaussNodes = 9
  nNodesPerElmQ1 = 4
  nNodesPerElmQ2 = 9

  w,s,t = GaussQuadPoints2D(3)

  tempArr4 = Array{Float64}(undef,nNodesPerElmQ1)
  tempArr9 = Array{Float64}(undef,nNodesPerElmQ2)

  xN     = copy(tempArr9); yN     = copy(tempArr9)

  phi    = copy(tempArr9); dphids = copy(tempArr9); dphidt=copy(tempArr9)
  dphidx = copy(tempArr9); dphidy = copy(tempArr9)

  psi    = copy(tempArr4); dpsids = copy(tempArr4); dpsidt = copy(tempArr4)
  dpsidx = copy(tempArr4); dpsidy = copy(tempArr4)

  # generate temporary element block matrices
  At = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Bt = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Ct = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ1)
  Dt = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Et = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Ft = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ1)  
  Gt = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)
  Ht = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)

  IA = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  JA = copy(IA); SA = copy(IA)
  IB = copy(IA); JB = copy(IA); SB = copy(IA)
  ID = copy(IA); JD = copy(IA); SD = copy(IA)
  IE = copy(IA); JE = copy(IA); SE = copy(IA)

  IG = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JG = copy(IG); SG = copy(IG)
  IH = copy(IG); JH = copy(IG); SH = copy(IG)

  # transpose!!!!
  IC = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JC = copy(IC); SC = copy(IC)
  IF = copy(IC); JF = copy(IC); SF = copy(IC)

  BigMat = spzeros(2*nUNodes+nPNodes,2*nUNodes+nPNodes)

  for el = 1:nElm
    localmat(mesh,el,xN,yN,w,s,t,
             At,Bt,Ct,Dt,Et,Ft,Gt,Ht,
             phi,dphidx,dphidy,dphids,dphidt,
             psi,dpsidx,dpsidy,dpsids,dpsidt,
             parameter...)
    # big block matrices -- A, B, D, E
    for ti = 1:nNodesPerElmQ2, tj = 1:nNodesPerElmQ2
      # At
      iii = mesh.cm[el].NodeList[ti]
      jjj = mesh.cm[el].NodeList[tj]
      BigMat[iii,jjj] += At[ti,tj]

      # Bt 
      iii = mesh.cm[el].NodeList[ti]
      jjj = mesh.cm[el].NodeList[tj] + nUNodes
      BigMat[iii,jjj] += Bt[ti,tj]   
      
      # Dt
      iii = mesh.cm[el].NodeList[ti] + nUNodes
      jjj = mesh.cm[el].NodeList[tj]
      BigMat[iii,jjj] += Dt[ti,tj]     

      # Et
      iii = mesh.cm[el].NodeList[ti] + nUNodes
      jjj = mesh.cm[el].NodeList[tj] + nUNodes
      BigMat[iii,jjj] += Et[ti,tj]  

      # At
      #IA[el,ti,tj] = mesh.cm[el].NodeList[ti]
      #JA[el,ti,tj] = mesh.cm[el].NodeList[tj]
      #SA[el,ti,tj] = At[ti,tj]
      # Bt
      #IB[el,ti,tj] = mesh.cm[el].NodeList[ti]
      #JB[el,ti,tj] = mesh.cm[el].NodeList[tj] + nUNodes
      #SB[el,ti,tj] = Bt[ti,tj]
      # Dt
      #ID[el,ti,tj] = mesh.cm[el].NodeList[ti] + nUNodes
      #JD[el,ti,tj] = mesh.cm[el].NodeList[tj]
      #SD[el,ti,tj] = Dt[ti,tj]
      # Et
      #IE[el,ti,tj] = mesh.cm[el].NodeList[ti] + nUNodes
      #JE[el,ti,tj] = mesh.cm[el].NodeList[tj] + nUNodes
      #SE[el,ti,tj] = Et[ti,tj]
    end
    # bottom block matrices -- G, H
    for ti = 1:nNodesPerElmQ1, tj = 1:nNodesPerElmQ2
      # Gt
      iii = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      jjj = mesh.cm[el].NodeList[tj]
      BigMat[iii,jjj] += Gt[ti,tj]

      # Ht
      iii = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      jjj = mesh.cm[el].NodeList[tj]  + nUNodes
      BigMat[iii,jjj] += Ht[ti,tj]

      # Ct
      iii = mesh.cm[el].NodeList[tj]
      jjj = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      BigMat[iii,jjj] += Ct[tj,ti]

      # Ft
      iii = mesh.cm[el].NodeList[tj]  + nUNodes
      jjj = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      BigMat[iii,jjj] += Ft[tj,ti]

      # Gt
      #IG[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      #JG[el,ti,tj] = mesh.cm[el].NodeList[tj]
      #SG[el,ti,tj] = Gt[ti,tj]
      # Ht
      #IH[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      #JH[el,ti,tj] = mesh.cm[el].NodeList[tj]  + nUNodes
      #SH[el,ti,tj] = Ht[ti,tj]
      # Ct
      #IC[el,ti,tj] = mesh.cm[el].NodeList[tj]
      #JC[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      #SC[el,ti,tj] = Ct[tj,ti]
      # Ft
      #IF[el,ti,tj] = mesh.cm[el].NodeList[tj]  + nUNodes
      #JF[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      #SF[el,ti,tj] = Ft[tj,ti]
    end
  end

  # bug could be in how this is vectorized...
  #I = vcat(vec(IA),vec(IB),vec(IC),vec(ID),vec(IE),vec(IF),vec(IG),vec(IH))
  #J = vcat(vec(JA),vec(JB),vec(JC),vec(JD),vec(JE),vec(JF),vec(JG),vec(JH))
  #S = vcat(vec(SA),vec(SB),vec(SC),vec(SD),vec(SE),vec(SF),vec(SG),vec(SH))

  #return sparse(I,J,S,2*nUNodes+nPNodes,2*nUNodes+nPNodes)
  return BigMat
end
"""
  assembleMPB(mesh,localmat,parameter...)

Assembles full 3x3 block matrix for fluid problems. Not sure how it differs 
  from the above.....
"""
function assembleMPB(mesh,localmat,parameter...)
  nElm = length(mesh.cm)
  nUNodes = length(mesh.xy); nPNodes = length(mesh.xyp)
  nGaussNodes = 9
  nNodesPerElmQ1 = 4
  nNodesPerElmQ2 = 9

  w,s,t = GaussQuadPoints2D(3)

  tempArr4 = Array{Float64}(undef,nNodesPerElmQ1)
  tempArr9 = Array{Float64}(undef,nNodesPerElmQ2)

  xN     = copy(tempArr9); yN     = copy(tempArr9)

  phi    = copy(tempArr9); dphids = copy(tempArr9); dphidt=copy(tempArr9)
  dphidx = copy(tempArr9); dphidy = copy(tempArr9)

  psi    = copy(tempArr4); dpsids = copy(tempArr4); dpsidt = copy(tempArr4)
  dpsidx = copy(tempArr4); dpsidy = copy(tempArr4)

  # generate temporary element block matrices
  At = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Ct = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ1)
  Et = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ2)
  Ft = Array{Float64}(undef,nNodesPerElmQ2,nNodesPerElmQ1)  
  Gt = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)
  Ht = Array{Float64}(undef,nNodesPerElmQ1,nNodesPerElmQ2)

  IA = zeros(Float64,nElm,nNodesPerElmQ2,nNodesPerElmQ2)
  JA = copy(IA); SA = copy(IA)
  IE = copy(IA); JE = copy(IA); SE = copy(IA)

  IG = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JG = copy(IG); SG = copy(IG)
  IH = copy(IG); JH = copy(IG); SH = copy(IG)

  # transpose!!!!
  IC = zeros(Float64,nElm,nNodesPerElmQ1,nNodesPerElmQ2)
  JC = copy(IC); SC = copy(IC)
  IF = copy(IC); JF = copy(IC); SF = copy(IC)

  for el = 1:nElm
    localmat(mesh,el,xN,yN,w,s,t,
             At,Ct,Et,Ft,Gt,Ht,
             phi,dphidx,dphidy,dphids,dphidt,
             psi,dpsidx,dpsidy,dpsids,dpsidt,
             parameter...)
    # big block matrices -- A, B, D, E
    for ti = 1:nNodesPerElmQ2, tj = 1:nNodesPerElmQ2
      # At
      IA[el,ti,tj] = mesh.cm[el].NodeList[ti]
      JA[el,ti,tj] = mesh.cm[el].NodeList[tj]
      SA[el,ti,tj] = At[ti,tj]
      # Et
      IE[el,ti,tj] = mesh.cm[el].NodeList[ti] + nUNodes
      JE[el,ti,tj] = mesh.cm[el].NodeList[tj] + nUNodes
      SE[el,ti,tj] = Et[ti,tj]
    end
    # bottom block matrices -- G, H
    for ti = 1:nNodesPerElmQ1, tj = 1:nNodesPerElmQ2
      # Gt
      IG[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      JG[el,ti,tj] = mesh.cm[el].NodeList[tj]
      SG[el,ti,tj] = Gt[ti,tj]
      # Ht
      IH[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      JH[el,ti,tj] = mesh.cm[el].NodeList[tj]  + nUNodes
      SH[el,ti,tj] = Ht[ti,tj]
      # Ct
      IC[el,ti,tj] = mesh.cm[el].NodeList[tj]
      JC[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      SC[el,ti,tj] = Ct[tj,ti]
      # Ft
      IF[el,ti,tj] = mesh.cm[el].NodeList[tj]  + nUNodes
      JF[el,ti,tj] = mesh.cmp[el].NodeList[ti] + 2*nUNodes
      SF[el,ti,tj] = Ft[tj,ti]
    end
  end

  # bug could be in how this is vectorized...
  I = vcat(vec(IA),vec(IC),vec(IE),vec(IF),vec(IG),vec(IH))
  J = vcat(vec(JA),vec(JC),vec(JE),vec(JF),vec(JG),vec(JH))
  S = vcat(vec(SA),vec(SC),vec(SE),vec(SF),vec(SG),vec(SH))

  return sparse(I,J,S,2*nUNodes+nPNodes,2*nUNodes+nPNodes)
end

"""
WeakScalar2D(mesh,farr)

Computes the weak form of some function farr using cartesian integral

ASSUMPTIONS:
- mesh is either Mesh or FluidMesh type

INPUT: 
  mesh:       [AbstractMesh]
                  mesh type that PDE is being solved on
  farr:       [Vector{Real}]
                  forcing values at mesh nodes

OUTPUT:
  F           [Vector{Real}]
                  RHS vector in form of RHS of linear Finite Element system
"""
function WeakScalar2D(mesh,farr)
  Nnodes     = length(mesh.xy)
  Nel        = length(mesh.cm)
  NumElNodes = length(mesh.cm[1].NodeList)
  F          = zeros(Nnodes)

  gaussorder = 0
  NumGaussNodes = 0
  if mesh.order == :Linear
    gaussorder = 2
    NumGaussNodes = Int(gaussorder^2)
  elseif mesh.order == :Quadratic
    gaussorder = 3
    NumGaussNodes = Int(gaussorder^2)
  end

  xNodes = zeros(Float64,NumGaussNodes)
  yNodes = zeros(Float64,NumGaussNodes)

  # Generate local F by elements
  # generate basis values at gauss points
  w,s,t = GaussQuadPoints2D(gaussorder)
  for el=1:Nel
    c = mesh.cm[el].NodeList
    for i=1:NumGaussNodes
      xNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].x
      yNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].y
    end

    for gpt=1:NumElNodes
      phi,_,_,jac = derivShape2D(s[gpt],t[gpt],xNodes,yNodes,gaussorder-1)

	    for i=1:NumElNodes
	      fcoefs      = farr[c]
	      fgpt        = shapeEval(fcoefs,phi)
	      F[mesh.cm[el].NodeList[i]] += fgpt*phi[i]*w[gpt]*jac
	    end
	  end
  end

  return F
end

"""
  WeakScalarAS(mesh,farr)

Computes the weak form of some function farr using axisymmetric integral

ASSUMPTIONS:
- mesh is either Mesh or FluidMesh type
- axisymmetric coordinates are accounted for. 
  Expects farr without modification by e.g. radial coordinate

INPUT: 
  mesh:       [AbstractMesh]
                  mesh type that PDE is being solved on
  farr:       [Vector{Real}]
                  forcing values at mesh nodes

OUTPUT:
  F           [Vector{Real}]
                  RHS vector in form of RHS of linear Finite Element system
                    in axisymmetric coordinates
"""
function WeakScalarAS(mesh,farr)
  Nnodes     = length(mesh.xy)
  Nel        = length(mesh.cm)
  NumElNodes = length(mesh.cm[1].NodeList)
  F          = zeros(Nnodes)

  gaussorder = 0
  NumGaussNodes = 0
  if mesh.order == :Linear
    gaussorder = 2
    NumGaussNodes = Int(gaussorder^2)
  elseif mesh.order == :Quadratic
    gaussorder = 3
    NumGaussNodes = Int(gaussorder^2)
  end

  xNodes = zeros(Float64,NumGaussNodes)
  yNodes = zeros(Float64,NumGaussNodes)

  # Generate local F by elements
  # generate basis values at gauss points
  w,s,t = GaussQuadPoints2D(gaussorder)
  for el=1:Nel
    c = mesh.cm[el].NodeList
    for i=1:NumGaussNodes
      xNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].x
      yNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].y
    end

    for gpt=1:NumElNodes
      phi,_,_,jac = derivShape2D(s[gpt],t[gpt],xNodes,yNodes,gaussorder-1)

	    for i=1:NumElNodes
	      fcoefs = farr[c]
        fgpt   = shapeEval(fcoefs,phi)
        yg     = shapeEval(yNodes,phi)
	      F[mesh.cm[el].NodeList[i]] += fgpt*phi[i]*w[gpt]*jac*yg
	    end
	  end
  end

  return F
end

# generates SUPG adjustment for axisymmetric cases
# parameter is of type for advection-diffusion, includes velocity field
# expects parameter to be a constant
function WeakSUPGAS(mesh,farr,parameter::T) where T<:AbstractConstantParameter
  Nnodes     = length(mesh.xy)
  Nel        = length(mesh.cm)
  NumElNodes = length(mesh.cm[1].NodeList)
  F          = zeros(Nnodes)

  # define 'wind'
  U = parameter.u
  V = parameter.v

  gaussorder = 0
  NumGaussNodes = 0
  if mesh.order == :Linear
    gaussorder = 2
    NumGaussNodes = Int(gaussorder^2)
  elseif mesh.order == :Quadratic
    gaussorder = 3
    NumGaussNodes = Int(gaussorder^2)
  end

  xNodes  = zeros(Float64,NumGaussNodes)
  yNodes  = zeros(Float64,NumGaussNodes)
  UNodes  = zeros(Float64,NumGaussNodes)
  VNodes  = zeros(Float64,NumGaussNodes)
  phi     = Array{Float64}(undef,NumElNodes)
  dphidx  = Array{Float64}(undef,NumElNodes)
  dphidy  = Array{Float64}(undef,NumElNodes)
  dphids  = Array{Float64}(undef,NumElNodes)
  dphidt  = Array{Float64}(undef,NumElNodes)


  # Generate local F by elements
  # generate basis values at gauss points
  w,s,t = GaussQuadPoints2D(gaussorder)
  for el=1:Nel
    c = mesh.cm[el].NodeList
    for i=1:NumGaussNodes
      xNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].x
      yNodes[i] = mesh.xy[mesh.cm[el].NodeList[i]].y
      UNodes[i] = U[c[i]]
      VNodes[i] = V[c[i]]
    end

    # compute h (max element length)
    hk = sqrt(areaCalc(xNodes,yNodes))

    # compute element Peclet number
    wmag = 0.0
    elPe = 0.0
    if mesh.order == :Linear
      uavg = sum(U[c])/4
      vavg = sum(V[c])/4
      wmag = sqrt(uavg^2 + vavg^2)
  
      elPe = 0.5*wmag*hk/parameter.κ
    elseif mesh.order == :Quadratic
      wmag = sqrt(U[c[9]]^2 + V[c[9]]^2)
  
      elPe = 0.5*wmag*hk/parameter.κ
    else 
      error("mesh.order not specified correctly")
    end

    # compute delta for SUPG
    δ = 0.5*hk/wmag*(1.0 - 1.0/elPe)

    for gpt=1:NumElNodes
      shape2D!(phi,dphids,dphidt,s[gpt],t[gpt],gaussorder-1)
      derivShape2D!(phi,dphidx,dphidy,dphids,dphidt,s[gpt],t[gpt],xNodes,yNodes,gaussorder-1)
      ug  = shapeEval(UNodes,phi)
      vg  = shapeEval(VNodes,phi)

	    for i=1:NumElNodes
	      fcoefs = farr[c]
        fgpt   = shapeEval(fcoefs,phi)
        yg     = shapeEval(yNodes,phi)
	      F[mesh.cm[el].NodeList[i]] += δ*fgpt*(ug*dphidx[i] + vg*dphidy[i])*w[gpt]*yg
	    end
	  end
  end

  return F
end

# given 4 nodes, computes area
function areaCalc(xNodes,yNodes)
  area1 = xNodes[1]*(yNodes[2] - yNodes[3]) + 
          xNodes[2]*(yNodes[3] - yNodes[1]) + 
          xNodes[3]*(yNodes[1] - yNodes[2])
  area2 = xNodes[1]*(yNodes[4] - yNodes[3]) + 
          xNodes[4]*(yNodes[3] - yNodes[1]) + 
          xNodes[3]*(yNodes[1] - yNodes[4])

  return 0.5*(abs(area1) + abs(area2))
end