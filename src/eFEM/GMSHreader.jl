# Code to be loaded into eFEMpart

abstract type AbstractFileType end

struct GMSHfile <: AbstractFileType
  file_name::String
end

# Define constructor for Mesh that takes in GMSH file
function Mesh(file_obj::GMSHfile)
  xy,cm,bndNodes,bndNames,order = GMSH_QUAD_READER(file_obj.file_name)

  bndNodes2 = GMSH_EXTRACTBOUNDARY_ALL(bndNodes,length(bndNames))
  bndNames2 = convert(Vector{Symbol},bndNames)
  bndNames3 = Index(bndNames2)
  bdry = DictList{Int}(bndNodes2,bndNames3)

  # generate element list
  cmEl = [Element(cm[i,:]) for i in 1:size(cm,1)]

  # generate node data structure from xy,xyp
  xyNodes = [Node(xy[i,1],xy[i,2]) for i=1:size(xy,1)]

  # turn order into a symbol (fast comparisons)
  order==1 ? o = Symbol("Linear") : o = Symbol("Quadratic")

  return ScalarMesh(xyNodes,cmEl,bdry,o)
end

# Define constructor for Fluid Mesh that takes in GMSH file
function FluidMesh(file_obj::GMSHfile)
  # generate normal mesh routine -- assumes that order of elements is 2
  xy,cm,bndNodes,bndNames,order = GMSH_QUAD_READER(file_obj.file_name)
  xyp,cmp = pMesh(xy,cm)

  bndNodes2 = GMSH_EXTRACTBOUNDARY_ALL(bndNodes,length(bndNames))
  bndNames2 = convert(Vector{Symbol},bndNames)
  bndNames3 = Index(bndNames2)
  bdry = DictList{Int}(bndNodes2,bndNames3)

  # generate element list
  cmEl  = [Element(cm[i,:])  for i in 1:size(cm,1)]
  cmpEl = [Element(cmp[i,:]) for i in 1:size(cmp,1)]

  # generate node data structure from xy,xyp
  xyNodes  = [Node(xy[i,1],xy[i,2])   for i=1:size(xy,1)]
  xypNodes = [Node(xyp[i,1],xyp[i,2]) for i=1:size(xyp,1)]

  # turn order into a symbol (fast comparisons)
  order==1 ? o = Symbol("Linear") : o = Symbol("Quadratic")

  return FluidMesh(xyNodes,cmEl,bdry,o,xypNodes,cmpEl)
end

"""
  printBlankLine(f)

prints blank line to IOstream handle f
"""
function printBlankLine(f)
  println(f,"")
end

"""
  GMSH_QUAD_READER(file_name) -> xy,cm,bndNodes

reading in GMSH file that assumes 2D quadrilateral elements.

To determine boundaries, you need to look at *.msh file header to see identifying numbers of boundaries, and afterwards group them how you want.
"""
function GMSH_QUAD_READER(file_name::String)
  # ---
  # open file
  # ---

  f          = open(file_name,"r")
	lines      = readlines(f)           # can be large ... might want to find
                                      # alternative method extract file info
	lineLength = length(lines)

	# ---
  # initialize tracking values
  # ---

  surfaceObjID = -1
  lineObjID    = -1
  npe          = -1

  numNodes = 0
  nodesLineStart = 0
  numElements = 0
  elementsLineStart = 0
  numBoundaries = 0
  numLineObjs = 0

  order = 0
  bndNames::Vector{String} = []

  # determine order
  for lineNumber = 1:lineLength
    ln = lines[lineNumber]
    if isElements(ln)
      elmStr = split(lines[lineNumber+2]," ")
      elType = str2int(elmStr[2])
      elType == 1 ? order=1 : order=2
    end
  end

  # ---
  # correct numbering expected for order
  # ---

  if order==1
    surfaceObjID = 3
    lineObjID = 1
    npe = 4
  elseif order==2
    surfaceObjID = 10
    lineObjID = 8
    npe = 9
  else
    #error("invalid order entered"::String)
  end

  # ---
  # First pass through file to obtain parameters
	# ---
	for lineNumber = 1:lineLength
	  ln = lines[lineNumber]
	  # boundary data
    if isBoundary(ln)
      numBoundaries = parse(Int,lines[lineNumber+1])

      for i=(lineNumber+2):(lineNumber+2+numBoundaries-2)
        str = lines[i]
        name = str[6:end-1]
        push!(bndNames,name)
      end
    end
    # node data
    if isNodes(ln)
      nodesLineStart = lineNumber+2
      numNodes = parse(Int,lines[lineNumber+1])
    end
    # element data
    if isElements(ln)
      elementsLineStart = lineNumber+2
      numElements = parse(Int,lines[lineNumber+1])

      for elm = 1:numElements
        elmStr = split(lines[lineNumber+elm+1]," ")
        elType = str2int(elmStr[2])
        if elType == lineObjID
          numLineObjs += 1
        end
      end
    end
  end

  # check that 'interior' surface exists
  ln = lines[5 + numBoundaries]
  name = ln[6:end-1]
  if !(isInterior(name))
    error("invalid mesh: requires PhysicalBoundary surface 'interior'")
  end

  # adjust element indexing to remove line objects
  numElements       -= numLineObjs
  lineObjsStart      = elementsLineStart
  elementsLineStart += numLineObjs

  # ---
  # generate arrays using above parameters
  # ---

  xy = zeros(Float64,(numNodes,2))
  cm = zeros(Int,(numElements,npe))
  # needs check in case numNodes < numElements, which happens for small meshes
  if numNodes < (numElements + numLineObjs)
    bndNodes = zeros(Int,(numNodes+numElements+numLineObjs,numBoundaries))
  else
    bndNodes = zeros(Int,(numNodes,numBoundaries))
  end

  # ---
  # obtain node data
  # ---
  nodeCounter = 1
  for lineNumber = nodesLineStart:(nodesLineStart+numNodes-1)
    ln = lines[lineNumber]
    POINTS = split(ln," ")
    xy[nodeCounter,1] = str2float(POINTS[2])
    xy[nodeCounter,2] = str2float(POINTS[3])
    nodeCounter += 1
  end

  # ---
  # obtain boundary data
  # ---
  bndCounter = 1
  for lineNumber = lineObjsStart:(lineObjsStart+numLineObjs-1)
    ln = lines[lineNumber]
    elmStr = split(ln," ")

    skip = parse(Int,elmStr[3]) + 3

    # boundary array
    physicalEntity = parse(Int,elmStr[4])
    bndNodes[bndCounter  ,physicalEntity] = str2int(elmStr[skip+1])
    bndNodes[bndCounter+1,physicalEntity] = str2int(elmStr[skip+2])
    if lineObjID==8   # quadratic elements have 3 line nodes
      bndNodes[bndCounter+2,physicalEntity] = str2int(elmStr[skip+3])
      bndCounter += 1
    end
    bndCounter += 2
  end

  # ---
  # obtain element connectivity data
  # ---
  elmCounter = 1

  for lineNumber = elementsLineStart:(elementsLineStart+numElements-1)
    ln = lines[lineNumber]
    elmStr = split(ln," ")

    skip = parse(Int,elmStr[3]) + 3

    for k=1:npe
      cm[elmCounter,k] = str2float(elmStr[skip+k])
    end
    elmCounter += 1
  end

  return xy,cm,bndNodes,bndNames,order
end

"""
  isBoundary(str)

Checks whether string is '\$PhysicalNames' function to be used when reading in *.msh files from gmsh. Used in GSMH_QUAD_READER
"""
function isBoundary(str)
    if VERSION==v"0.5.0"
      str = str[1:end-1]
    end
    if str == "\$PhysicalNames"
        return true
    else
        return false;
    end
end

"""
  isElements(str)

Checks whether string is '\$Elements' function to be used when reading in *.msh files from gmsh. Used in GSMH_QUAD_READER
"""
function isElements(str)
    if VERSION==v"0.5.0"
      str = str[1:end-1]
    end
    if str == "\$Elements"
        return true
    else
        return false;
    end
end

"""
  isNodes(str)

Checks whether string is '\$Nodes' function to be used when reading in
*.msh files from gmsh. Used in GSMH_QUAD_READER
"""
function isNodes(str)
    if VERSION==v"0.5.0"
      str = str[1:end-1]
    end
    if str == "\$Nodes"
        return true
    else
        return false;
    end
end

function isInterior(str)
  if VERSION==v"0.5.0"
    str = str[1:end-1]
  end
  if str == "interior"
      return true
  else
      return false;
  end
end

"""
  str2int(str)

parses String into Int
"""
function str2int(str)
  return parse(Int,str)
end

"""
  str2float(str)

parses String into Float64
"""
function str2float(str)
  return parse(Float64,str)
end

"""
  GMSH_BOUNDARYPARSE(N)

Takes in gmsh boundary matrix and returns first N boundaries (sorted and remove 0 element).

Really crappy programming, just a stopgap until a more elegant solution is found
"""
function GMSH_EXTRACTBOUNDARY(bndNodes,i)
  # weird aspect of GSMH reader implementation includes 0 at beginning of sorted array. Need to include check in case implementation is changed later on
  return sort(unique(bndNodes[:,i]))[2:end]
end

"""
  GMSH_EXTRACTBOUNDARY_ALL(bndNodes,total)

completely stupid way to extract boundary node lists, simply used to shorten code length in driver files
"""
function GMSH_EXTRACTBOUNDARY_ALL(bndNodes,total)
  # weird aspect of GSMH reader implementation includes 0 at beginning of sorted array. Need to include check in case implementation is changed later on
  if total==1
    A = sort(unique(bndNodes[:,1]))[2:end]
    return [A]
  elseif total==2
    A = sort(unique(bndNodes[:,1]))[2:end]
    B = sort(unique(bndNodes[:,2]))[2:end]
    return [A,B]
  elseif total==3
    A = sort(unique(bndNodes[:,1]))[2:end]
    B = sort(unique(bndNodes[:,2]))[2:end]
    C = sort(unique(bndNodes[:,3]))[2:end]
    return [A,B,C]
  elseif total==4
    A = sort(unique(bndNodes[:,1]))[2:end]
    B = sort(unique(bndNodes[:,2]))[2:end]
    C = sort(unique(bndNodes[:,3]))[2:end]
    D = sort(unique(bndNodes[:,4]))[2:end]
    return [A,B,C,D]
  elseif total==5
    A = sort(unique(bndNodes[:,1]))[2:end]
    B = sort(unique(bndNodes[:,2]))[2:end]
    C = sort(unique(bndNodes[:,3]))[2:end]
    D = sort(unique(bndNodes[:,4]))[2:end]
    E = sort(unique(bndNodes[:,5]))[2:end]
    return [A,B,C,D,E]
  else
    error("total=$(total) is too large for this file")
  end
end
