# File to be loaded into eFEMpart

"""
  DomainNorm(xy,cm,u;normID='Inf')
  
Computes either the L-1, L-2, or L-inf norm on a given domain. Expects quadrilateral elements, and order 1 or 2 polynomial basis.

ASSUMPTIONS:
- Uses strings for L-norm specification ... change this to symbols

INPUT:
  xy:     [scalar array]
            Node List for node positions
  cm:     [Matrix]
            Connectivity Matrix
  u:      [Array]
            Scalar field to put in integrand
  normID: [Symbol]   
            indicates L-1,2,Inf norm
OUTPUT:
  norm:  [scalar array]
            (1,2,Inf)-Norm of field u
"""
function DomainNorm(xy,cm,u;normID=:Inf)
  norm  = 0.0
  Nel   = length(cm)
  ngpt  = length(cm[1].NodeList)
  order = (ngpt==4 ? 1 : 2)
  
  # generate basis values at gauss points
  w,s,t = GaussQuadPoints2D(order+1)
    
  # L1 or L2 norm
  if (normID==:1 || normID==:2)
    normIDnum = str2float(normID)
    
    for el=1:Nel
      xNodes = [xy[cm[el].NodeList[i]].x for i=1:ngpt]
      yNodes = [xy[cm[el].NodeList[i]].y for i=1:ngpt]
      
      for gpt=1:ngpt
        phi,_,_,jac = derivShape2D(s[gpt],t[gpt],xNodes,yNodes,order)
	      
	      ucoefs = [u[cm[el].NodeList[i]] for i=1:ngpt]
	      ugpt = shapeEval(ucoefs,phi)
	      
	      norm += w[gpt]*abs(ugpt)^normIDnum*jac
      end
    end

    # take sqrt of norm value for L2 norm    
    if normID==:2
      norm = sqrt(norm)
    end
    
  # infinity norm (aka max norm)
  elseif normID==:Inf
    norm = maximum(abs,u)
  else
    error("invalid normID given. Choices are '1', '2', or 'Inf'.")    
  end
  
  return norm
end 


"""
  hCalc(mesh)
  
Computes the spacing parameter h for a given mesh. Works on quadrilateral meshes

ASSUMPTIONS:
- Most effective when all elements are roughly the same size. This
    algorithm will not measure coarseness of stretched meshes accurately.

INPUT:
  mesh:  [AbstractMesh]
            Mesh that includes nodes and connectivity

OUTPUT:
  h:     [AbstractFloat]
            Maximum edge distance of mesh, measure for mesh refinement
"""
function hCalc(mesh)
  Nel = length(mesh.cm)
  h   = 0.0
  
  for el=1:Nel
    for i=1:4
      j = mod(i,4)+1  
      ci = mesh.cm[el].NodeList[i]
      cj = mesh.cm[el].NodeList[j]
      
      xdiff = mesh.xy[ci].x - mesh.xy[cj].x
      ydiff = mesh.xy[ci].y - mesh.xy[cj].y 
       
      htemp = sqrt(xdiff^2 + ydiff^2)  
      
      if htemp > h
        h = htemp
      end
            
    end
  end
  
  return h
end
