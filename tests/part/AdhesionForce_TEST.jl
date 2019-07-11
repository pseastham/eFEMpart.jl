# tests AdhesionForce functions

# Functions for computing particle-wall contact forces of 
# via LJ potentials and a quadrature

# loads in Point2D type
include("../../src/part/CellLists.jl")
# load wall types
include("../../src/part/ParticleTypes.jl")
# loads in LJ force function
include("../../src/part/CohesionForce.jl")
# load onSegment
include("../../src/part/isInside.jl")

# linewall
function NearestPoint!_linewall_TEST()
    point = Point2D(0.0,0.0)
    node = Point2D(0.0,1.0)
    linestart = Point2D(-1.0,0.0)
    lineend   = Point2D(1.0,0.0)
    wall = LineWall([linestart,lineend])

    #@code_warntype NearestPoint!(point,node,wall)
    @btime NearestPoint!($point,$node,$wall)
    

    nothing
end
# circlewall
function NearestPoint!_circlewall_TEST()
    point = Point2D(0.0,0.0)
    node = Point2D(0.0,1.0)
    center = Point2D(1.4,0.0)
    radius = 0.2
    wall = CircleWall(center,radius)

    #@code_warntype NearestPoint!(point,node,wall)
    @btime NearestPoint!($point,$node,$wall)
    
    println(point)

    nothing
end
# arcwall
function NearestPoint!_arcwall_TEST()
    point  = Point2D(0.0,0.0)
    node   = Point2D(0.0,1.0)

    r1     = Point2D(1.6,0.0)
    center = Point2D(1.4,0.0)
    r2     = Point2D(1.4,0.2)
    wall   = ArcWall(r1,center,r2)

    #@code_warntype NearestPoint!(point,node,wall)
    @btime NearestPoint!($point,$node,$wall)

    nothing
end

function isCloseEnough_TEST()
    node  = Point2D(0.0,1.0)
    point = Point2D(0.0,0.0)
    radius = 0.1
    k = 5.0
    rm = 1.5

    #@code_warntype isCloseEnough(point,node,k,rm,radius)
    @btime isCloseEnough($point,$node,$k,$rm,$radius)

    nothing
end

# linewall test
function generateQuadNodes!_linewall_TEST()
    node = Point2D(0.0,1.0)
    linestart = Point2D(-1.0,0.0)
    lineend   = Point2D(1.0,0.0)
    wall = LineWall([linestart,lineend])
    point = NearestPoint(node,wall)

    Nquad = 20
    k = 5.0
    rm = 1.5

    # initialize arrays
    xquad = zeros(Nquad); yquad = zeros(Nquad)

    #@code_warntype generateQuadNodes!(xquad,yquad,point,wall,node,k,rm)
    @btime generateQuadNodes!($xquad,$yquad,$point,$wall,$k,$rm)

    nothing
end
# circlewall
function generateQuadNodes!_circlewall_TEST()
    node   = Point2D(0.0,1.0)
    center = Point2D(1.4,0.0)
    radius = 0.2
    wall = CircleWall(center,radius)
    point = NearestPoint(node,wall)

    Nquad = 20
    k = 5.0
    rm = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    #@code_warntype generateQuadNodes!(xquad,yquad,point,wall,k,rm)
    @btime generateQuadNodes!($xquad,$yquad,$point,$wall,$k,$rm)

    nothing
end
# arcwall
function generateQuadNodes!_arcwall_TEST()
    node   = Point2D(0.0,1.0)
    r1     = Point2D(1.6,0.0)
    center = Point2D(1.4,0.0)
    r2     = Point2D(1.4,0.2)
    wall   = ArcWall(r1,center,r2)
    point = NearestPoint(node,wall)

    Nquad = 20
    k = 5.0
    rm = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    #@code_warntype generateQuadNodes!(xquad,yquad,point,wall,k,rm)
    @btime generateQuadNodes!($xquad,$yquad,$point,$wall,$k,$rm)

    nothing
end

# linewall
function isInLine_linewall_TEST()
    node = Point2D(0.0,1.0)
    linestart = Point2D(-1.0,0.0)
    lineend   = Point2D(1.0,0.0)
    wall = LineWall([linestart,lineend])
    point = NearestPoint(node,wall)

    Nquad = 20
    k = 5.0
    rm = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    generateQuadNodes!(xquad,yquad,point,wall,k,rm)

    # initialize vectors so they don't have to keep being re-allocated
    s = Point2D(0.0,0.0)

    #@code_warntype isInLine(wall,xquad[2],yquad[2],s)
    @btime isInLine($wall,$xquad[2],$yquad[2],$s) 

    nothing
end
# circlewall
function isInLine_circlewall_TEST()
    node   = Point2D(0.0,1.0)
    center = Point2D(1.4,0.0)
    radius = 0.2
    wall = CircleWall(center,radius)
    point = NearestPoint(node,wall)

    Nquad = 20
    k = 5.0
    rm = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    generateQuadNodes!(xquad,yquad,point,wall,k,rm)

    # initialize vectors so they don't have to keep being re-allocated
    s = Point2D(0.0,0.0)

    #@code_warntype isInLine(wall,xquad[2],yquad[2],s)
    @btime isInLine($wall,$xquad[2],$yquad[2],s) 

    nothing
end
# arcwall
function isInLine_arcwall_TEST()
    node   = Point2D(0.0,1.0)
    r1     = Point2D(1.6,0.0)
    center = Point2D(1.4,0.0)
    r2     = Point2D(1.4,0.2)
    wall   = ArcWall(r1,center,r2)
    point = NearestPoint(node,wall)

    Nquad = 20
    k = 5.0
    rm = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    generateQuadNodes!(xquad,yquad,point,wall,k,rm)

    # initialize vectors so they don't have to keep being re-allocated
    s = Point2D(0.0,0.0)

    #@code_warntype isInLine(wall,xquad[2],yquad[2],s)
    @btime isInLine($wall,$xquad[2],$yquad[2],$s) 
    
    nothing
end

# linewall
function wallTrapQuad_linewall_TEST()
    node = Point2D(0.0,0.2)
    linestart = Point2D(-1.0,0.0)
    lineend   = Point2D(1.0,0.0)
    wall = LineWall([linestart,lineend])
    point = NearestPoint(node,wall)

    particleradius = 0.1
    Nquad = 20
    k = 5.0
    rm = 1.5            # whats difference between rm and rc???
    ϵ = 1.0
    rc = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    generateQuadNodes!(xquad,yquad,point,wall,k,rm)

    s = Point2D(0.0,0.0)

    fx,fy = wallTrapQuad(node,wall,xquad,yquad,ϵ,particleradius,rc,s)
    #@code_warntype wallTrapQuad(node,wall,xquad,yquad,ϵ,particleradius,rc,s)
    @btime wallTrapQuad($node,$wall,$xquad,$yquad,$ϵ,$particleradius,$rc,$s)
    
    return fx,fy
end

# circle wall
function wallTrapQuad_circlewall_TEST()
    node   = Point2D(1.1,0.0)
    center = Point2D(1.4,0.0)
    radius = 0.2
    wall = CircleWall(center,radius)
    point = NearestPoint(node,wall)

    particleradius = 0.1
    Nquad = 20
    k = 5.0
    rm = 1.5            # whats difference between rm and rc???
    ϵ = 1.0
    rc = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    generateQuadNodes!(xquad,yquad,point,wall,k,rm)

    # initialize vectors so they don't have to keep being re-allocated
    s = Point2D(0.0,0.0)

    fx,fy = wallTrapQuad(node,wall,xquad,yquad,ϵ,particleradius,rc,s)
    #@code_warntype wallTrapQuad(node,wall,xquad,yquad,ϵ,particleradius,rc,s)
    @btime wallTrapQuad($node,$wall,$xquad,$yquad,$ϵ,$particleradius,$rc,$s)
    
    return fx,fy
end

# arc wall
function wallTrapQuad_arcwall_TEST()
    node   = Point2D(1.4,0.0)
    r1     = Point2D(1.6,0.0)
    center = Point2D(1.4,0.0)
    r2     = Point2D(1.4,0.2)
    wall   = ArcWall(r1,center,r2)
    point = NearestPoint(node,wall)

    particleradius = 0.1
    Nquad = 20
    k = 5.0
    rm = 1.5            # whats difference between rm and rc???
    ϵ = 1.0
    rc = 1.5

    # initialize vectors to be assigned in-place
    xquad = zeros(Nquad)
    yquad = zeros(Nquad)

    generateQuadNodes!(xquad,yquad,point,wall,k,rm)

    # initialize vectors so they don't have to keep being re-allocated
    s = Point2D(0.0,0.0)

    fx,fy = wallTrapQuad(node,wall,xquad,yquad,ϵ,particleradius,rc,s)
    #@code_warntype wallTrapQuad(node,wall,xquad,yquad,ϵ,particleradius,rc,s)
    @btime wallTrapQuad($node,$wall,$xquad,$yquad,$ϵ,$particleradius,$rc,$s)
    
    return fx,fy
end