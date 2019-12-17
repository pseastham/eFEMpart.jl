using eFEMpart
include("MeshTransform.jl")

point = [rand(),rand()]

p1 = 2*[-1.0,-1.0]
p2 = 2*[1.0,-1.0]
p3 = 2*[1.0,1.0]
p4 = 2*[-1.0,1.0]
quad = [p1,p2,p3,p4]

s,t = reverseQuadMap(quad,point)

