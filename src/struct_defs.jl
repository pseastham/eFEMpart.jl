mutable struct scratch_data
    n_particles::Int
    polygon::Vector{Point2D{Float64}}
    extremePoint::Point2D{Float64}
    pointOnWall::Point2D{Float64}
    w::Vector{Float64}
    uEl::Vector{Float64}
    vEl::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
    c::Vector{Float64}
    d::Vector{Float64}
    cfX::Vector{Float64}
    cfY::Vector{Float64}
    afX::Vector{Float64}
    afY::Vector{Float64}
    gfX::Vector{Float64}
    gfY::Vector{Float64}
    sfX::Vector{Float64}
    sfY::Vector{Float64}

    function scratch_data(n_particles::Int)
        polygon      = [Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0)]
        extremePoint = Point2D(100_000.0,0.0)
        pointOnWall  = Point2D(0.0,0.0)
        w   = zeros(4)
        uEl = zeros(4)
        vEl = zeros(4)
        a   = zeros(3)
        b   = zeros(3)
        c   = zeros(3)
        d   = zeros(3)
        cfX = zeros(n_particles)
        cfY = zeros(n_particles)
        afX = zeros(n_particles)
        afY = zeros(n_particles)
        gfX = zeros(n_particles)
        gfY = zeros(n_particles)
        sfX = zeros(n_particles)
        sfY = zeros(n_particles)
        return new(n_particles,polygon,extremePoint,pointOnWall,
                    w,uEl,vEl,a,b,c,d,
                    cfX,cfY,afX,afY,gfX,gfY,sfX,sfY)
    end
end  
