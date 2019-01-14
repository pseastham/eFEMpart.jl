# Code to be loaded into eFEMpart

# This file contains functions to convert the discrete
# sand particles into a continuous porosity field discretized onto
# a FEM mesh

# Finds area of quadrilateral using Heron's formula for triangle areas
function QuadArea(quad::Vector{Vector{Float64}})
    a = sideLength(quad[1],quad[2])
    b = sideLength(quad[2],quad[3])
    c = sideLength(quad[3],quad[1])
    A1 = HeronTriangleArea(a,b,c)
    a = sideLength(quad[1],quad[3])
    b = sideLength(quad[3],quad[4])
    c = sideLength(quad[4],quad[1])
    A2 = HeronTriangleArea(a,b,c)

    Area = A1+A2

    return Area
end

function QuadArea(mesh,el::Int)
    c = mesh.cm[el].NodeList[1:4]
    p1 = [mesh.xy[c[1]].x,mesh.xy[c[1]].y]
    p2 = [mesh.xy[c[2]].x,mesh.xy[c[2]].y]
    p3 = [mesh.xy[c[3]].x,mesh.xy[c[3]].y]
    p4 = [mesh.xy[c[4]].x,mesh.xy[c[4]].y]
    quad  = [p1,p2,p3,p4]

    return QuadArea(quad)
end

# computes area of triangle with length a,b,c
# using Heron's formula
function HeronTriangleArea(a::T,b::T,c::T) where T<:Real
    p = (a+b+c)/2
    A = sqrt(p*(p-a)*(p-b)*(p-c))
    return A
end

# find length of side given 2 points, p1 and p2
function sideLength(p1::T,p2::T) where T<:Real
    sl = sqrt((p2[1]-p1[1])^2 + (p2[2]-p1[2])^2)
    return sl
end

# computes areas of all elements in mesh
function ElementAreas(mesh::AbstractMesh)
    Nelms = length(mesh.cm)
    elArea  = zeros(Float64,Nelms)
    for el=1:Nelms
        elArea[el] = QuadArea(mesh,el::Int)
    end

    return elArea
end

function SandToPorosityGaussian(mesh,femCL,pList,rm,σ,wList)
  Nnodes     = length(mesh.xy)
  Nparticles = length(pList)
  femPList = [Particle(mesh.xy[i].x,mesh.xy[i].y) for i=1:Nnodes]

  porosity = zeros(Float64,Nnodes)

  for j=1:Nparticles
    SandToPorosityGaussianSingleParticle!(pList[j],mesh,rm,femCL,σ,porosity,wList)
  end

  return porosity
end

function SandToPorosityGaussianSingleParticle!(particle,mesh,rm,femCL,σ,porosity,wList)
    H = 100.0
    Ncells = length(femCL)
    pcc = 0             # cell index containing particle
    ϵ  = 1.0
    AreaOfParticle = pi*(0.5*rm)^2

    # find cell containing particle
    for i=1:Ncells
        rect = femCL[i].square
        point = [particle.xpos,particle.ypos]
        if isInsideRect(rect,point)
            pcc = i
            break
        end
    end

    # makes sure particle is within some cell
    if pcc != 0
        # array of neighbors & itself particle cell
        nbs = vcat(pcc,femCL[pcc].NeighborList)

        # find particles in all relevant cells
        pclose = Array{Int,1}(0)
        for i=1:length(nbs)
            ncl = length(femCL[nbs[i]].ParticleList)
            if ncl > 0
                for j=1:ncl
                    push!(pclose,femCL[nbs[i]].ParticleList[j])
                end
            end
        end

        # compute gaussian distance and circle area approx for particle
        for j in pclose
            p1 = [particle.xpos,particle.ypos]
            q1 = [mesh.xy[j].x,mesh.xy[j].y]
            # check that path doens't intersect wall
            if noIntersectingWall(p1,q1,wList)
                Δx = p1[1] - q1[1]
                Δy = p1[2] - q1[2] #particle.ypos - mesh.xy[j].y
                r = sqrt(Δx^2 + Δy^2)
                porosity[j] += GaussEval(r,σ,H)*AreaOfParticle
            end
        end
    end

    return nothing
end

# evaluates a 1D gaussian
function GaussEval(r::T,σ::T,H::T) where T<:Real
    ex = exp(-0.5*(r/σ)^2)
    return H*ex
end
