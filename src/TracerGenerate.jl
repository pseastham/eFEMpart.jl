# to be loaded into eFEM
# code to generate tracers for flow fields

abstract type AbstractTracerInfo end
abstract type AbstractTracerSource end

mutable struct Tracer
    xpos::Float64
    ypos::Float64
end

mutable struct TracerInfo 
    Δt::Float64                     # time step
    T::Float64                      # final time
    nT::Int                         # number of time steps
    TracerArr::Vector{Tracer}
end

struct TracerLineSource <: AbstractTracerSource
    LineStart::Vector{Float64}
    LineEnd::Vector{Float64}
    nTracers::Int
end

struct TracerBoxSource <: AbstractTracerSource

end

function GenerateTracers(mesh::FluidMesh,u,v,TI::TracerInfo,file_name)
    NTracers = length(TI.TracerArr)

    # Generate time array
    timeArr = 0.0:TI.Δt:TI.T

    # needed to check whether point is inside domain
    xyMesh = NodesToArray(mesh.xy)
    cmMesh = ElementsToArray(mesh.cm)
    
    xy = zeros(NTracers,2)
    removeArrayVals = []

    # initialize xy array
    for tj=1:NTracers
        xy[tj,1] = TI.TracerArr[tj].xpos
        xy[tj,2] = TI.TracerArr[tj].ypos
    end

    # export vtk save file
    PARTICLES_TO_VTK(xy,1.0,string(file_name,"_0");time=0.0)

    # for each time pointer
    for ti=1:TI.nT
        for tj=1:NTracers
            # check whether any nodes fell outside of domain on last entry
            if !(isInsideDomain(xyMesh,cmMesh,[TI.TracerArr[tj].xpos,TI.TracerArr[tj].ypos]))
                push!(removeArrayVals,tj)
            end
        end

        # check to generate new tracers -- not yet

        # interpolate velocity at tracer positions
        for tj=1:NTracers
            # checks that tracer is still in domain
            if !(tj in removeArrayVals)
                # move tracers by flow field (requires interpolation)
                uval = pointTransform(mesh,u,[TI.TracerArr[tj].xpos,TI.TracerArr[tj].ypos])
                vval = pointTransform(mesh,v,[TI.TracerArr[tj].xpos,TI.TracerArr[tj].ypos])

                TI.TracerArr[tj].xpos += uval*TI.Δt
                TI.TracerArr[tj].ypos += vval*TI.Δt

                xy[tj,1] = TI.TracerArr[tj].xpos
                xy[tj,2] = TI.TracerArr[tj].ypos
            end
        end

        # export vtk save file
        PARTICLES_TO_VTK(xy,1.0,string(file_name,"_",ti);time=ti*TI.Δt)
    end
    nothing    
end