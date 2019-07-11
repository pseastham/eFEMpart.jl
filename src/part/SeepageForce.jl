# define and test barycentric velocity interpolation, to be used
# for computing seepage force in particle simulations

include("isInside.jl")
include("CellLists.jl")
include("BarycentricInterpolation.jl")

"""
computes velocity of single particle without using cell list

efficient in that it allocates no new memory, but inefficient in that
it loops over all elements (no cell list given)
"""
function BarycentricVelocityInterp(mesh,u::Vector{T},v::Vector{T},particle::Point2D,
                                    polygon::Vector{Point2D},extremeArr::Point2D,
                                    w::Vector{T},uEl::Vector{T},vEl::Vector{T},
                                    a::Vector{T},b::Vector{T},c::Vector{T},
                                    d::Vector{T}) where T<:Real
    sVx = 0.0; sVy = 0.0

    elInd = 0
    extremeArr.y = particle.y

    # find which element the particle belongs to
    for el=1:length(mesh.cm)
        # construct polygon
        for ti=1:4
            polygon[ti].x = mesh.xy[mesh.cm[el].NodeList[ti]].x
            polygon[ti].y = mesh.xy[mesh.cm[el].NodeList[ti]].y
        end

        if isInside(polygon,4,particle;extreme=extremeArr)
            elInd = el
            break
        end
    end

    # find weights for barycentric  
    computeBaryWeights!(w,polygon,particle,a,b,c,d)

    # interpolate values 
    for ti=1:4
        uEl[ti] = u[mesh.cm[elInd].NodeList[ti]]
        vEl[ti] = v[mesh.cm[elInd].NodeList[ti]]
    end

    sVx = dot(uEl, w)
    sVy = dot(vEl, w)

    return sVx,sVy
end

"""
efficient velocity interpolation.

Given nodelist of particles (via point positions), interpolates velocity from FEM using 
cell list.

for max efficiency, velocities are interpolated according to cell ordering, NOT particle
index ordering

passes in a vector for velocity interpolations (expected u and v are same size!!)
"""
function BarycentricVelocityInterp_CL!(uInterp::Vector{Float64},vInterp::Vector{Float64},mesh,u::Vector{T},v::Vector{T},
                                        nodeList::Vector{Point2D},polygon::Vector{Point2D},particleCL::CellList,
                                        meshCLmap::Array{Array{Int64,N} where N,1},extremePoint::Point2D,w::Vector{T},
                                        uEl::Vector{T},vEl::Vector{T},a::Vector{T},
                                        b::Vector{T},c::Vector{T},d::Vector{T}) where T<:Real
    if !(length(uInterp) == length(vInterp) == length(nodeList))
        throw(DimensionMismatch("uInterp, vInterp and nodeList need to be all same length"))
    end

    # loop over cells
    for cellInd = 1:length(particleCL.cells)
        for nodeInd in particleCL.cells[cellInd].nodeList
            foundNode = false

            extremePoint.y = nodeList[nodeInd].y

            # check whether cell belongs to element INSIDE cell
            for elInd in meshCLmap[cellInd]
                # construct polygon -- FEM element
                for ti=1:4
                    polygon[ti].x = mesh.xy[mesh.cm[elInd].NodeList[ti]].x
                    polygon[ti].y = mesh.xy[mesh.cm[elInd].NodeList[ti]].y
                end

                if isInside(polygon,4,nodeList[nodeInd];extreme=extremePoint)
                    foundNode = true
                    # find weights for barycentric  
                    computeBaryWeights!(w,polygon,nodeList[nodeInd],a,b,c,d)

                    # interpolate values 
                    for ti=1:4
                        uEl[ti] = u[mesh.cm[elInd].NodeList[ti]]
                        vEl[ti] = v[mesh.cm[elInd].NodeList[ti]]
                    end

                    uInterp[nodeInd] = dot(uEl, w)
                    vInterp[nodeInd] = dot(vEl, w)

                    break
                end
            end
            
            # last resort in case particle wasn't found using cell list
            if !(foundNode)
                uInterp[nodeInd],vInterp[nodeInd] = BarycentricVelocityInterp(mesh,u,v,nodeList[nodeInd],polygon,extremePoint,
                                                                                w,uEl,vEl,a,b,c,d)
                @warn "particle found in cell list -- must use brute-force velocity interpolation" maxlog=1
            end
        end
    end

    nothing
end