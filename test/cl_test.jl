using eFEMpart

# test fem_generate_map
rect = [-1.0,1.0,-1.0,1.0]
N    = 64
mesh = squareMesh(rect,N,2)
L    = 0.1
totalBounds = [rect[1]-L,rect[2]+L,rect[3]-L,rect[4]+L]
femCLmap = StokesParticles.femGenerateMap(mesh,totalBounds,L)

# test getElementFromNode!
rect = [-1.0,1.0,-1.0,1.0]
N    = 5
mesh = squareMesh(rect,N,1)
elementArray = zeros(Int,4)
nodeInd = 12
@test StokesParticles.getElementFromNode!(elementArray,mesh,nodeInd) != π     # just to make sure function runs