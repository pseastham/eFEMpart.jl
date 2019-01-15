# File to be loaded into eFEMpart

function meshgrid(vx::AbstractVector{T}, vy::AbstractVector{T}) where T
	m, n = length(vy), length(vx)
  vx = reshape(vx, 1, n)
  vy = reshape(vy, m, 1)

  return (repmat(vx, m, 1)', repmat(vy, 1, n)')
end

# square [x0, x1, y0, y1]
function squareMesh(square,N,order)
  x = linspace(square[1],square[2],N)
  y = linspace(square[3],square[4],N)

  # create nodes
  X,Y = meshgrid(x,y)
  xy = [X[:] Y[:]]

  # create element connectivity
  ni = size(X,1); # number of rows
  nNodes = size(xy,1);
  t2nidxMap = 1:nNodes-ni;
  topNode = ni:ni:nNodes-ni;
  t2nidxMap = deleteat!(collect(t2nidxMap),collect(topNode))
  k = t2nidxMap;
  cm = [k k+1 k+ni+1 k+ni]

  if order==2
    xy,cm = Order1ToOrder2(xy,cm)
    xy,cm = deleteDuplicates(xy,cm)
  end

  # create boundary nodes
  left=Int[]; top=Int[]; right=Int[]; bottom=Int[]
  for i=1:size(xy,1)
    # left
    if isMachZero(xy[i,1],square[1])
      push!(left,i)
    end
    # right
    if isMachZero(xy[i,1],square[2])
      push!(right,i)
    end
    # bottom
    if isMachZero(xy[i,2],square[3])
      push!(bottom,i)
    end
    # top
    if isMachZero(xy[i,2],square[4])
      push!(top,i)
    end
  end

  # generate bdry object
  bndNames = Index([:left,:top,:right,:bottom])
  bdry = DictList{Int}([left,top,right,bottom],bndNames)

  (order==1 ? o=:Linear : o=:Quadratic)

  # Generate Mesh!
  return ScalarMesh(ArrayToNodes(xy),ArrayToElements(cm),bdry,o)
end

function squareMeshFluid(square,N)
  mesh = squareMesh(square,N,2)

  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  xyptemp,cmptemp = pMesh(xy,cm)

  xyp = ArrayToNodes(xyptemp)
  cmp = ArrayToElements(cmptemp)

  return FluidMesh(mesh.xy,mesh.cm,mesh.bdry,:Quadratic,xyp,cmp)
end

function axisymMesh(Nr,Nt,Rmax,order;phi=(x -> x))
  xi = linspace(0,1,Nr)
  t = linspace(0,pi,Nt)

  # create nodes
  XI,THETA = meshgrid(xi,t)
  xy = [XI[:] THETA[:]]

  # create element connectivity
  ni = size(XI,1); # number of rows
  nNodes = size(xy,1)
  t2nidxMap = 1:nNodes-ni;
  topNode = ni:ni:nNodes-ni;
  t2nidxMap = deleteat!(collect(t2nidxMap),collect(topNode))
  k = t2nidxMap;
  cm = [k k+1 k+ni+1 k+ni]

  if order==2
    xy,cm = Order1ToOrder2(xy,cm)
    xy,cm = deleteDuplicates(xy,cm)
  end

  # create boundary nodes
  inflow=Int[]; outflow=Int[]; axis=Int[]; sphere=Int[]
  for i=1:size(xy,1)
    # axis
    if isMachZero(xy[i,2],0.0) || isMachZero(xy[i,2],Float64(pi))
      push!(axis,i)
    end
    # inflow
    if isMachZero(xy[i,1],1) && xy[i,2]<=Float64(pi)/2
      push!(inflow,i)
    end
    # outflow
    if isMachZero(xy[i,1],1) && xy[i,2]>Float64(pi)/2
      push!(outflow,i)
    end
    # sphere
    if isMachZero(xy[i,1],0.0)
      push!(sphere,i)
    end
  end

  # transform xy nodes
  transform(r) = exp(log(Rmax)*phi(r))

  XI2    = xy[:,1]
  THETA2 = xy[:,2]

  R = transform.(XI2)
  X = R.*cos.(THETA2)
  Y = R.*sin.(THETA2)
  xy = [X[:] Y[:]]

  # generate bdry object
  bndNames = Index([:inflow,:outflow,:axis,:sphere])
  bdry = DictList{Int}([inflow,outflow,axis,sphere],bndNames)

  (order==1 ? o=:Linear : o=:Quadratic)

  # Generate Mesh!
  return ScalarMesh(ArrayToNodes(xy),ArrayToElements(cm),bdry,o)
end

function axisymFluid(Nr,Nt,Rmax;phi=(x -> x))
  mesh = axisymMesh(Nr,Nt,Rmax,2;phi=phi)

  xy = NodesToArray(mesh.xy)
  cm = ElementsToArray(mesh.cm)

  xyptemp,cmptemp = pMesh(xy,cm)

  xyp = ArrayToNodes(xyptemp)
  cmp = ArrayToElements(cmptemp)

  return FluidMesh(mesh.xy,mesh.cm,mesh.bdry,:Quadratic,xyp,cmp)
end

# function to convert array to Nodes
function ArrayToNodes(xy)
  return [Node(xy[i,1],xy[i,2]) for i=1:size(xy,1)]
end

# function to convert Nodes to array form
function NodesToArray(xy)
    # get xy into correct form
  x = [xy[i].x for i=1:length(xy)]
  y = [xy[i].y for i=1:length(xy)]
  return hcat(x,y)
end

# function to convert array to elements
function ArrayToElements(cm)
  return [Element(cm[i,:]) for i in 1:size(cm,1)]
end

# function to convert Elements to array form
function ElementsToArray(cm)
  return [cm[i].NodeList[j] for i=1:length(cm), j=1:length(cm[1].NodeList)]
end

# converts order 1 meshes to order 2 meshes
function Order1ToOrder2(xy,cm)
  nNodesOld = size(xy,1)
  nElm      = size(cm,1)

  xNew = Float64[]
  yNew = Float64[]
  cmNew = Array{Int}(nElm,9)

  cmNew[:,1:4] = cm   # will throw error if cm isn't quad

  # for each element
  newNodeInd = nNodesOld + 1
  for el=1:nElm
    # for nodes 5 - 9
    for i=5:9
      x=0.0; y=0.0
      if i==5
        x = (xy[cm[el,1],1] + xy[cm[el,2],1])/2
        y = (xy[cm[el,1],2] + xy[cm[el,2],2])/2
      elseif i==6
        x = (xy[cm[el,2],1] + xy[cm[el,3],1])/2
        y = (xy[cm[el,2],2] + xy[cm[el,3],2])/2
      elseif i==7
        x = (xy[cm[el,3],1] + xy[cm[el,4],1])/2
        y = (xy[cm[el,3],2] + xy[cm[el,4],2])/2
      elseif i==8
        x = (xy[cm[el,4],1] + xy[cm[el,1],1])/2
        y = (xy[cm[el,4],2] + xy[cm[el,1],2])/2
      elseif i==9
        x = (xy[cm[el,1],1] + xy[cm[el,2],1] +
             xy[cm[el,3],1] + xy[cm[el,4],1])/4
        y = (xy[cm[el,1],2] + xy[cm[el,2],2] +
             xy[cm[el,3],2] + xy[cm[el,4],2])/4
      end
      # add node to xy
      push!(xNew,x)
      push!(yNew,y)

      # add node index to cm
      cmNew[el,i] = newNodeInd
      newNodeInd += 1
    end
  end

  xyNew = Array{Float64}(newNodeInd-1,2)

  for i=1:(newNodeInd-1)
    if i <= nNodesOld
      xyNew[i,1] = xy[i,1]; xyNew[i,2] = xy[i,2]
    else
      k = i-nNodesOld
      xyNew[i,1] = xNew[k]; xyNew[i,2] = yNew[k]
    end
  end

  return xyNew, cmNew
end

# need to make function that replaces/deletes duplicate nodes
# for the love of all that is holy, do not change anything in here. the indices
# are a nightmare
function deleteDuplicates(xy,cm)
  nNodesOld = size(xy,1)
  nElm   = size(cm,1)

  nodeIndPos5 = view(cm,:,5)
  nodeIndPos6 = view(cm,:,6)
  nodeIndPos7 = view(cm,:,7)
  nodeIndPos8 = view(cm,:,8)

  xyNew = Float64[]
  skipList = Int[]
  cmNew = copy(cm)

  # for each node in 5 or 7 position
  for ti=1:length(nodeIndPos5)
    i = nodeIndPos5[ti]
    if !(i in skipList)
      for tj=1:length(nodeIndPos7)
        j = nodeIndPos7[tj]
        # check whether x AND y match within epsilon
        if isMachZero(xy[i,1],xy[j,1]) &&
           isMachZero(xy[i,2],xy[j,2]) &&
           i!=j
          # find other index for node j
          kjfind = findfirst(cm,j)
          # 1) assign corresponding node index to cm of other node
          cmNew[kjfind] = i
          # 3) add other node to list so that it isn't crossed later on
          push!(skipList,j)
          break
        end
      end
    end
  end

  # for each node in 6 or 8 position
  for ti=1:length(nodeIndPos6)
    i = nodeIndPos6[ti]
    if !(i in skipList)
      for tj=1:length(nodeIndPos8)
        j = nodeIndPos8[tj]
        # check whether x AND y match within epsilon
        if isMachZero(xy[i,1],xy[j,1]) &&
           isMachZero(xy[i,2],xy[j,2]) &&
           i!=j
          # find other index for node j
          kjfind = findfirst(cmNew,j)
          # 1) assign corresponding node index to cm of other node
          cmNew[kjfind] = i
          # 3) add other node to list so that it isn't crossed later on
          push!(skipList,j)
          break
        end
      end
    end
  end

  cmsubtract = zeros(Int,(length(cm)))
  # forward sweep to decide how many to subtract from each cm
  tempind=1
  for i in cmNew
    for k in skipList
      if i>=k
        cmsubtract[tempind]+=1
      end
    end
    tempind+=1
  end

  # backward sweep to actually take off those deleted indices
  for i=1:length(cm)
    cmNew[i] -= cmsubtract[i]
  end

  xyNew = xy[setdiff(1:end, skipList), :]

  return xyNew, cmNew
end

# function to check whether two numbers are within machine epsilon
function isMachZero(x,y)
  ep = eps()
  if abs(x-y) < ep
    return true
  else
    return false
  end
end

"""
  pMesh(xy,cm) -> xyp, cmp

Generates pressure mesh equivalents given a velocity mesh (xy,cm).

Might need to generate map between velocity and pressure nodes ....?
"""
function pMesh(xy,cm)
  Nel         = size(cm,1)
  NumPNodes   = length(unique(@view cm[:,1:4]))
  xyp         = zeros(NumPNodes,2)
  cmp         = zeros(Int,(Nel,4))

  nodeCheckarr::Array{Int,1} = []

  elCounter = 1
  nodeCounter = 1
  nodeCheck = 1
  for el=1:Nel
    for i=1:4
      c = cm[el,i]
      if c in nodeCheckarr
        cmp[elCounter,i] = findfirst(nodeCheckarr,c)
      else
        push!(nodeCheckarr,c)
        xyp[nodeCounter,:] = @view xy[c,:]
        cmp[elCounter,i] = nodeCounter
        nodeCounter += 1
      end

      nodeCheck += 1
    end

    elCounter += 1
  end

  return xyp,cmp
end


"""
  lowerUV(xy,cm,u,v) -> u,v

computes velocity on coarser mesh using p-mesh
"""
function lowerUV(xy,cm,u,v)
  Nel       = size(cm,1)
  NumPNodes = length(unique(cm[:,1:4]))
  q1u   = zeros(NumPNodes,1)
  q1v   = zeros(NumPNodes,1)
  cmp = 0*Array{Int}(Nel,4)
  nodeCheckarr::Array{Int,1} = []

  elCounter = 1
  nodeCounter = 1
  nodeCheck = 1
  for el=1:Nel
    for i=1:4
      c = cm[el,i]
      if c in nodeCheckarr
        cmp[elCounter,i] = findfirst(nodeCheckarr,c)
      else
        push!(nodeCheckarr,c)
        q1u[nodeCounter] = u[c]
        q1v[nodeCounter] = v[c]
        nodeCounter += 1
      end

      nodeCheck += 1
    end

    elCounter += 1
  end

  return q1u,q1v
end

"""
  lowerScalar(xy,cm,q2s) -> q1s

computes scalar on coarser mesh using p-mesh
"""
function lowerScalar(xy,cm,q2s)
  Nel       = size(cm,1)
  NumQ1Nodes = length(unique(cm[:,1:4]))
  q1s   = zeros(NumQ1Nodes,1)
  #cmp = 0*Array{Int}(Nel,4)
  nodeCheckarr::Array{Int,1} = []

  elCounter = 1
  nodeCounter = 1
  nodeCheck = 1
  for el=1:Nel
    for i=1:4
      c = cm[el,i]
      if c in nodeCheckarr
        #cmp[elCounter,i] = findfirst(nodeCheckarr,c)
      else
        push!(nodeCheckarr,c)
        q1s[nodeCounter] = q2s[c]
        nodeCounter += 1
      end

      nodeCheck += 1
    end

    elCounter += 1
  end

  return q1s
end
