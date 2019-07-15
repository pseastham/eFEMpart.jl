# define and test barycentric velocity interpolation

using BenchmarkTools: @btime

include("../../src/part/SeepageForce.jl")  # loads functions to be tested

function interp_timing(N=5)
    # generate mesh
    Dom = [-1.0,1.0,-1.0,1.0]
    #N    = 5
    mesh = squareMesh(Dom,N,2)

    # generate u,v data
    u = zeros(length(mesh.xy))
    v = zeros(length(mesh.xy))
    for ti=1:length(mesh.xy)
        u[ti] = 1.0 - mesh.xy[ti].y^2 
    end    

    # generate particle position
    xpos=0.999; ypos=0.999
    p = Point2D(xpos,ypos)

    # print out FEM stats
    #FEMstats(mesh)

    # initialize memory so that allocation isn't an issue
    polygon    = [Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0)]
    extremeArr = Point2D(100_000.0,0.0)
    w          = zeros(4)
    uEl        = zeros(4)
    vEl        = zeros(4)
    a = zeros(3)
    b = zeros(3)
    c = zeros(3)
    d = zeros(3)

    # call 
    #@code_warntype BarycentricVelocityInterp(mesh,u,v,p,polygon,extremeArr)
    #@btime BarycentricVelocityInterp($mesh,$u,$v,$p,$polygon,$extremeArr,$w,$uEl,$vEl,$a,$b,$c,$d)
    sVx,sVy = BarycentricVelocityInterp(mesh,u,v,p,polygon,extremeArr,w,uEl,vEl,a,b,c,d)

    if true
        println()
        println("expected u:    ",1-ypos^2)
        println("inteprolate u: ",sVx)
        println()
        println("expected v:    ",0.0)
        println("inteprolate v: ",sVy)
        println()
    end

    nothing
end

# timing of everything the same but using a cell list
function interp_timing_withCL(N=5,L=0.1)
    # generate mesh
    Dom = [-1.0,1.0,-1.0,1.0]
    #N    = 5
    mesh = squareMesh(Dom,N,2)

    # generate u,v data
    u = zeros(length(mesh.xy))
    v = zeros(length(mesh.xy))
    for ti=1:length(mesh.xy)
        u[ti] = 1.0 - mesh.xy[ti].y^2 
    end

    # generate particle position
    xpos=1.0; ypos=1.0
    p = Point2D(xpos,ypos)

    # generate node cell list with p
    nodeList = [p]
    totalBounds = [-1.0,1.0,-1.0,1.0]
    particleCL = generateCellList(nodeList,totalBounds,L)
    
    # generate mesh cell list map
    meshCLmap = femGenerateMap(mesh,totalBounds,L)

    # initialize memory so that allocation isn't an issue
    polygon    = [Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0),Point2D(0.0,0.0)]
    extremeArr = Point2D(100_000.0,0.0)
    w          = zeros(4)
    uEl        = zeros(4)
    vEl        = zeros(4)
    a          = zeros(3)
    b          = zeros(3)
    c          = zeros(3)
    d          = zeros(3)
    uInterp    = zeros(1)
    vInterp    = zeros(1)

    # call 
    #@code_warntype BarycentricVelocityInterp_CL!(mesh,u,v,nodeList,polygon,particleCL,meshCLmap,extremeArr,w,uEl,vEl,a,b,c,d,uInterp,vInterp)
    @btime BarycentricVelocityInterp_CL!($uInterp,$vInterp,$mesh,$u,$v,$nodeList,$polygon,$particleCL,$meshCLmap,$extremeArr,$w,$uEl,$vEl,$a,$b,$c,$d)
    BarycentricVelocityInterp_CL!(uInterp,vInterp,mesh,u,v,nodeList,polygon,particleCL,meshCLmap,extremeArr,w,uEl,vEl,a,b,c,d)

    if false
        println()
        println("expected u:    ",1-ypos^2)
        println("inteprolate u: ",uInterp[1])
        println()
        println("expected v:    ",0.0)
        println("inteprolate v: ",vInterp[1])
        println()
    end

    # testing
    return mesh
end