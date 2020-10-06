using Parameters

@with_kw struct Stokes_Para
    # mesh
    Nintervals::Int = 32        # number of intervals for your mesh
    xmin::Float64 = -2.0        # minimum x-value for mesh
    xmax::Float64 = 2.0         # minimum x-value for mesh
    ymin::Float64 = -2.0        # minimum x-value for mesh
    ymax::Float64 = 1.0         # minimum x-value for mesh 
    
    # parameters
    μ::Float64 = 1.0            # viscosity

    # boundary conditions
    lBCval::Float64 = 1.0       # strength of inflow

    # visualization
    solFolder::String = "solutions"         # folder for solution to go in
    solFile::String = "stokes_example"     # name of vtk solutions output
end

# test comment