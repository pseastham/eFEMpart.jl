# Defines new types to be used in particle simulations

# ==============================================================================
# TYPE DEFINITIONS
# ==============================================================================
abstract type AbstractWall end

mutable struct Point2D
    x::Float64
    y::Float64
end

mutable struct LineWall <: AbstractWall
    nodes::Vector{Point2D}        # 2 points, defining start and end
    n::Vector{Float64}            # normal unit vector of wall
    t::Vector{Float64}            # tangent unit vector of wall
    orientation::Symbol           # can be :right, :left, or :both to determine 
end                               #     on which side of the circle the domain
                                  #     of the problem is. Helpful when computing
                                  #     interacting forces. :right means to the 
                                  #     right if standing on start and looking at end. 
                                  #     :left is opposite, and :both is both.

mutable struct CircleWall <: AbstractWall
    center::Point2D               # point that defines center
    radius::Float64               # radius of circle
    orientation::Symbol           # can be :inward, :outward, or :both to determine 
end                               #     on which side of the circle the domain
                                  #     of the problem is. Helpful when computing
                                  #     interacting forces. :inward is towards center
                                  #     of circle

mutable struct ArcWall <: AbstractWall
    nodes::Vector{Point2D}        # 3 nodes, 1) start, 2) center, and 3) end going CCW (counterclockwise)
    orientation::Symbol           # can be :inward, :outward, or :both to determine 
end                               #     on which side of the circle the domain
                                  #     of the problem is. Helpful when computing
                                  #     interacting forces. :inward is towards center
                                  #     of circle

# ==============================================================================
# INITIALIZATION DEFINITIONS
# ==============================================================================

# initializes new wall defined by 2 points (Start,End)
function LineWall(nodes::Vector{Point2D},orientation::Symbol)
    if length(nodes) != 2
        throw(DimensionMismatch("line is defined by 2 points"))
    end

    x0 = [nodes[1].x,nodes[2].x]
    x1 = [nodes[1].y,nodes[2].y]

    n = zeros(Float64,2)
    t = zeros(Float64,2)

    # compute wall length
    WL = sqrt((nodes[2].x - nodes[1].x)^2 + (nodes[2].y - nodes[1].y)^2)

    # comute tangent
    t[1] = (nodes[2].x - nodes[1].x)/WL
    t[2] = (nodes[2].y - nodes[1].y)/WL

    n[1] = -t[2]
    n[2] =  t[1]

    return LineWall(nodes,n,t,orientation)
end
## allows for entering points by entering lines
#function LineWall(V::Vector{Vector{Float64}},orientation)
#    p1 = Point(V[1][1],V[1][2])
#    p2 = Point(V[2][1],V[2][2])
#    return LineWall([p1,p2],orientation)
#end
# default orientation is :both
LineWall(vORpoint) = LineWall(vORpoint,:both)

# Default orientation is :both
CircleWall(V,r) = CircleWall(V,r,:both)

function ArcWall(V::Vector{Vector{Float64}},orientation::Symbol)
    p1 = Point(V[1][1],V[1][2])
    p2 = Point(V[2][1],V[2][2])
    p3 = Point(V[3][1],V[3][2])

    r1 = sqrt((p2.x-p1.x)^2 + (p2.y-p1.y)^2)
    r2 = sqrt((p2.x-p3.x)^2 + (p2.y-p3.y)^2)

    TOL = 1e-12
    if abs(r1-r2) > TOL
        error("points do not make a circle with constant radius, r1=$(round(r1,2)), r2=$(round(r2,2))")
    end

    return ArcWall([p1,p2,p3],orientation)
end
# Default Orientation is :both
ArcWall(V) = ArcWall(V,:both)
ArcWall(v1,v2,v3) = ArcWall([v1,v2,v3])
