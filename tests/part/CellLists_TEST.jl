# Generates and updates cell list object, and generates mapping for 
# finite element mesh elements to a corresponding cell list cells

using BenchmarkTools

include("../../src/part/ParticleTypes.jl")  # to load in Point2D type

function isInsideRect_TEST()
  sampleRect = [0.0, 1.0, 0.0, 1.0]
  samplePt = Point2D(rand(),rand())

  @btime isInsideRect($sampleRect,$samplePt)
end

function generateCellList_TEST()
  N = 1_000
  nodeList = [Point2D(rand(),rand()) for i=1:N]
  TotalBounds = [0.0,1.0,0.0,1.0]
  L = 0.1

  #@code_warntype generateCellList(nodeList,TotalBounds,L)
  @btime generateCellList($nodeList,$TotalBounds,$L)

  cellList = generateCellList(nodeList,TotalBounds,L)

  sum = 0 
  for i=1:length(cellList.cells)
    sum += length(cellList.cells[i].nodeList)
  end
  println(N," = ",sum," ?")

  return cellList, nodeList
end

function updateCellList!_TEST()
  N = 1_000
  # generate cell list
  nodeList = [Point2D(rand(),rand()) for i=1:N]
  TotalBounds = [0.0,1.0,0.0,1.0]
  L = 0.1

  cl = generateCellList(nodeList,TotalBounds,L)

  # greate new cell list
  nodeList2 = [Point2D(rand(),rand()) for i=1:N]
  @code_warntype updateCellList!(cl,nodeList2)
  #@btime updateCellList!(deepcopy($(cl)),$nodeList2)

  sum = 0 
  for i=1:length(cl.cells)
    sum += length(cl.cells[i].nodeList)
  end
  println(sum)

  nothing
end

"""
  not actually faster
"""
function updateCellList_SMART!_TEST()
  N = 1_000
  # generate cell list
  nodeList = [Point2D(rand(),rand()) for i=1:N]
  TotalBounds = [0.0,1.0,0.0,1.0]
  L = 0.1

  cl = generateCellList(nodeList,TotalBounds,L)

  # greate new cell list
  for i=1:N
    nodeList[i].x += 2*L*(2*rand() - 1.0)
    nodeList[i].y += 2*L*(2*rand() - 1.0)
  end
  @btime updateCellList_SMART!(deepcopy($(cl)),$nodeList)

  sum = 0 
  for i=1:cl.nCells
    sum += length(cl.cells[i].nodeList)
  end
  println(sum)

  nothing
end

function femGenerateMap_TEST()
  # generate mesh
  Dom = [-1.0,1.0,-1.0,1.0]
  N    = 64
  mesh = squareMesh(Dom,N,2)

  L = 0.1
  totalBounds = [Dom[1]-L,Dom[2]+L,Dom[3]-L,Dom[4]+L]

  #@code_warntype FEMgenerateMap(mesh,totalBounds,L)
  @btime femGenerateMap($mesh,$totalBounds,$L)
  femCLmap = femGenerateMap(mesh,totalBounds,L)

  println(typeof(femCLmap))

  return femCLmap
end

function getElementFromNode!_TEST()
    # generate mesh
    Dom = [-1.0,1.0,-1.0,1.0]
    N    = 5
    mesh = squareMesh(Dom,N,1)

    elementArray = zeros(Int,4)

    nodeInd = 12

    #@code_warntype getElementFromNode!(elementArray,mesh,nodeInd)
    @btime getElementFromNode!($elementArray,$mesh,$nodeInd)
    #getElementFromNode!(elementArray,mesh,nodeInd)

    if false
      println("node #",nodeInd)
      println("elm array: ",elementArray)
      for ti=1:4
        if elementArray[ti] != 0
          println("cm",elementArray[ti],": ",mesh.cm[elementArray[ti]].NodeList)
        end
      end
    end

    nothing
end