# File to be loaded into SinkholeModule

# allows for perfect reflection when particle encounters a wall
# OBSOLETE with new way of governing interactions, namely, with LJ 
#   forces instead of contact forces
function WallBounce(wall::LineWall,vec::Vector{Float64})
    tx = wall.t[1]
    ty = wall.t[2]
    nx = wall.n[1]
    ny = wall.n[2]

    vt = vec[1]*tx + vec[2]*ty
    vn = vec[1]*nx + vec[2]*ny

    newvx = vt*wall.t[1] - vn*wall.n[1]
    newvy = vt*wall.t[2] - vn*wall.n[2]

    return newvx,newvy
end

# checks whether a particle is inside a wall or not
# OBSOLETE with new way of defining walls
function checkInsideWall(particle,WallList)
    Nwalls = length(WallList)

    for wallIndex=1:Nwalls
        polygon = [[WallList[wallIndex].brect[i].x,
                    WallList[wallIndex].brect[i].y] for i=1:4]
        if isInside(polygon, 4, [particle.xpos,particle.ypos])
            return true, wallIndex
        end
    end

    return false, 0
end

# checks whether particle-to-evaluation-point path (p1,q1) crosses wall (p2,q2)
function noIntersectingWall(p1,q1,WallList)
    Nwalls = length(WallList)

    for i=1:Nwalls
        p2 = [WallList[i].nodes[1].x,WallList[i].nodes[1].y]
        q2 = [WallList[i].nodes[2].x,WallList[i].nodes[2].y]
        if doIntersect(p1,q1,p2,q2)
            return false
        end
    end

    return true
end

