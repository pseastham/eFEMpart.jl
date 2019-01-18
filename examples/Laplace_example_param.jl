using Parameters

@with_kw struct Lap_Para
    # mesh
    Nintervals::Int = 32        # number of intervals for your mesh
    basisorder::Int = 2         # element basis function order (1 or 2)
    xmin::Float64 = -2.0        # minimum x-value for mesh
    xmax::Float64 = 2.0         # minimum x-value for mesh
    ymin::Float64 = -1.0        # minimum x-value for mesh
    ymax::Float64 = 1.0         # minimum x-value for mesh 
    
    # boundary conditions
    lBCval::Float64 = 5.0       # (dirichlet) value at left boundary
    rBCval::Float64 = 0.0       # (dirichlet) value at right boundary

    # visualization
    solFolder::String = "solutions"         # folder for solution to go in
    solFile::String = "laplace_example"     # name of vtk solutions output
    solName::String = "temperature"         # name of unknown
end