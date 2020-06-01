# define and test barycentric velocity interpolation, to be used
# for computing seepage force in particle simulations

using LinearAlgebra: dot

"""
computes velocity of single particle without using cell list

efficient in that it allocates no new memory, but inefficient in that
it loops over all elements (no cell list given)
"""
function compute_seepage_force_nocl!(sfX::Vector{T},sfY::Vector{T},mesh,u::Vector{T},v::Vector{T},pList,data) where T<:Real
    fill!(sfX,zero(T))
    fill!(sfY,zero(T))

    for ti=1:data.n_particles
        elFound = false
        elInd = 0
        data.extremePoint.y = pList[ti].pos.y

        # find which element the particle belongs to
        for el=1:length(mesh.cm)
            # construct polygon
            for ti=1:4
                data.polygon[ti].x = mesh.xy[mesh.cm[el].NodeList[ti]].x
                data.polygon[ti].y = mesh.xy[mesh.cm[el].NodeList[ti]].y
            end

            if StokesParticles.is_inside(data.polygon,4,pList[ti].pos;extreme=data.extremePoint)
                elInd = el
                elFound = true
                break
            end
        end

        # find weights for barycentric  
        computeBaryWeights!(data.w,data.polygon,pList[ti].pos,data.a,data.b,data.c,data.d)

        if elFound
            # interpolate values 
            for tj=1:4
                data.uEl[tj] = u[mesh.cm[elInd].NodeList[tj]]
                data.vEl[tj] = v[mesh.cm[elInd].NodeList[tj]]
            end

            sfX[ti] = dot(data.uEl, data.w)
            sfY[ti] = dot(data.vEl, data.w)
        else
            sfX[ti] = 0.0
            sfY[ti] = 0.0
        end
    end

    nothing
end

function compute_single_seepage_force_nocl(mesh,u::Vector{T},v::Vector{T},particle,data) where T<:Real
    sfX = 0.0
    sfY = 0.0
    elFound = false
    elInd = 0
    data.extremePoint.y = particle.pos.y

    # find which element the particle belongs to
    for el=1:length(mesh.cm)
        # construct polygon
        for ti=1:4
            data.polygon[ti].x = mesh.xy[mesh.cm[el].NodeList[ti]].x
            data.polygon[ti].y = mesh.xy[mesh.cm[el].NodeList[ti]].y
        end

        if StokesParticles.is_inside(data.polygon,4,particle.pos;extreme=data.extremePoint)
            elInd = el
            elFound = true
            break
        end
    end

    # find weights for barycentric  
    computeBaryWeights!(data.w,data.polygon,particle.pos,data.a,data.b,data.c,data.d)

    if elFound
        # interpolate values 
        for tj=1:4
            data.uEl[tj] = u[mesh.cm[elInd].NodeList[tj]]
            data.vEl[tj] = v[mesh.cm[elInd].NodeList[tj]]
        end

        sfX = dot(data.uEl, data.w)
        sfY = dot(data.vEl, data.w)
    end

    return sfX, sfY
end

"""
efficient velocity interpolation.

Given nodelist of particles (via point positions), interpolates velocity from FEM using 
cell list.

for max efficiency, velocities are interpolated according to cell ordering, NOT particle
index ordering

passes in a vector for velocity interpolations (expected u and v are same size!!)
"""
function compute_seepage_force!(sfX::Vector{T},sfY::Vector{T},mesh,u::Vector{T},v::Vector{T},pList,data,cl,femap) where T<:Real
    fill!(sfX,zero(T))
    fill!(sfY,zero(T))

    # loop over cells
    for cellInd = 1:length(cl.cells)
        for nodeInd in cl.cells[cellInd].particleIDList
            foundNode = false

            data.extremePoint.y = pList[nodeInd].pos.y

            # check whether cell belongs to element INSIDE cell
            for elInd in femap[cellInd]
                # construct polygon -- FEM element
                for ti=1:4
                    data.polygon[ti].x = mesh.xy[mesh.cm[elInd].NodeList[ti]].x
                    data.polygon[ti].y = mesh.xy[mesh.cm[elInd].NodeList[ti]].y
                end

                if StokesParticles.is_inside(data.polygon,4,pList[nodeInd].pos;extreme=data.extremePoint)
                    foundNode = true
                    # find weights for barycentric  
                    computeBaryWeights!(data.w,data.polygon,pList[nodeInd].pos,data.a,data.b,data.c,data.d)

                    # interpolate values 
                    for ti=1:4
                        data.uEl[ti] = u[mesh.cm[elInd].NodeList[ti]]
                        data.vEl[ti] = v[mesh.cm[elInd].NodeList[ti]]
                    end
                    data.sfX[nodeInd] = dot(data.uEl, data.w)
                    data.sfY[nodeInd] = dot(data.vEl, data.w)
                    break
                end
            end

            if !(foundNode)
                sfX[nodeInd],sfY[nodeInd] = compute_single_seepage_force_nocl(mesh,u,v,pList[nodeInd],data)
            end
        end
    end

    nothing
end