using Parameters

@with_kw struct AD_Para
    # mesh
    Nintervals::Int = 50        # number of intervals for your mesh
    basisorder::Int = 2         # element basis function order (1 or 2)
    xmin::Float64 = -1.0        # minimum x-value for mesh
    xmax::Float64 = 1.0         # minimum x-value for mesh
    ymin::Float64 = -1.0        # minimum x-value for mesh
    ymax::Float64 = 1.0         # minimum x-value for mesh 
    
    # boundary conditions
    dirVal::Float64 = 1.0       # (dirichlet) value at left boundary
    Pe::Float64 = 20.0           # Peclet number (advection / diffusion)

    # visualization
    solFolder::String = "solutions"         # folder for solution to go in
    solFile::String = "ad_example"     # name of vtk solutions output
    solName::String = "chemical"         # name of unknown
end