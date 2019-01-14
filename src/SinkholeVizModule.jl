# Module used to Plot various aspects of
# Kept in seperate module to avoid loading time of Plots during normal
#   run-time

module SinkholeVizModule

using eFEMpart, Plots, LaTeXStrings
font = Plots.font("DejaVu Sans", 14)
pyplot(label="",color="black",size=(800,800),border=true,
        guidefont=font, xtickfont=font, ytickfont=font, legendfont=font,
        markersize=8,linewidth=3,ratio=:equal)

export plotPoints,
       plotWall,
       plotWalls,
       plotParticles!

function plotWalls(wallList::Vector{T},rm::Float64) where T<:AbstractWall
  N = length(wallList)

  p = plot() #initialize plot
  for i=1:N
    plotWall(wallList[i],p,rm)
  end
  
  return p
end

function plotParticles!(plotObj,pList,rm)
    colorInd = [:black,:red,:green,:yellow,:purple,:orange,:blue]
    N = length(pList)

    for i=1:N
        tc = mod(i,length(colorInd))+1
        plotParticle!(pList[i],plotObj,colorInd[tc],rm)
    end
    
    return plotObj
end

function plotWall(wall::LineWall,plotObj::Plots.Plot,rm)
    radius = 0.5*rm

    x1 = wall.nodes[1].x; x2 = wall.nodes[2].x
    y1 = wall.nodes[1].y; y2 = wall.nodes[2].y
    θWALL = atan2(y2-y1,x2-x1)

    Nθ = 15
    θ1 = linspace(-pi/2 + θWALL,  pi/2 + θWALL,Nθ)
    θ2 = linspace( pi/2 + θWALL,3*pi/2 + θWALL,Nθ)

    len = sqrt((x1-x2)^2 + (y1-y2)^2)
    
    tx = (x2-x1)/len
    ty = (y2-y1)/len

    nx = -ty
    ny = tx

    # 1 arc
    x = [x2 + radius*cos(θ1[i]) for i=1:Nθ]
    y = [y2 + radius*sin(θ1[i]) for i=1:Nθ]
    plot!(plotObj,x,y)

    # 2 line
    x = [x2 + radius*nx,x1 + radius*nx]
    y = [y2 + radius*ny,y1 + radius*ny]
    plot!(plotObj,x,y)

    # 3 arc 
    x = [x1 + radius*cos(θ2[i]) for i=1:Nθ]
    y = [y1 + radius*sin(θ2[i]) for i=1:Nθ]
    plot!(plotObj,x,y)
    
    # 4 line
    x = [x1 - radius*nx,x2 - radius*nx]
    y = [y1 - radius*ny,y2 - radius*ny]
    plot!(plotObj,x,y)
end
function plotWall(wall::CircleWall,plotObj::Plots.Plot,rm)
    radius = 0.5*rm
    Nplot = 50; Δθ = 2*pi/(Nplot-1)
    curvex = [wall.center.x + (wall.radius + radius)*cos((i-1)*Δθ) for i=1:Nplot]
    curvey = [wall.center.y + (wall.radius + radius)*sin((i-1)*Δθ) for i=1:Nplot]
    plot!(plotObj,curvex,curvey)
end
function plotWall(wall::ArcWall,plotObj::Plots.Plot,rm)
    radius = 0.5*rm

    pa = [wall.nodes[1].x,wall.nodes[1].y]
    c  = [wall.nodes[2].x,wall.nodes[2].y]
    pb = [wall.nodes[3].x,wall.nodes[3].y]
    r = sqrt((pa[1]-c[1])^2 + (pa[2]-c[2])^2)
    θ1 = atan2(pa[2]-c[2],pa[1]-c[1])
    θ2 = atan2(pb[2]-c[2],pb[1]-c[1])

    if θ1 > θ2
        θ2 += 2*pi
    end

    Nplot = 50; Δθ = (θ2 - θ1)/(Nplot-1)

    # 1 big arc -- done
    curvex = [c[1] + (r+radius)*cos((i-1)*Δθ+θ1) for i=1:Nplot]
    curvey = [c[2] + (r+radius)*sin((i-1)*Δθ+θ1) for i=1:Nplot]
    plot!(plotObj,curvex,curvey)

    # 2 far cap
    θWALL = θ2 + pi/2
    θ2cap = linspace(-pi/2 + θWALL,  pi/2 + θWALL,Nplot)
    x = [pb[1] + radius*cos(θ2cap[i]) for i=1:Nplot]
    y = [pb[2] + radius*sin(θ2cap[i]) for i=1:Nplot]
    plot!(plotObj,x,y)

    # 3 small arc 
    curvex = [c[1] + (r-radius)*cos((i-1)*Δθ+θ1) for i=1:Nplot]
    curvey = [c[2] + (r-radius)*sin((i-1)*Δθ+θ1) for i=1:Nplot]
    plot!(plotObj,curvex,curvey)
    
    # 4 close cap
    θWALL = θ1 + pi/2
    θ1cap = linspace( pi/2 + θWALL,3*pi/2 + θWALL,Nplot)
    x = [pa[1] + radius*cos(θ1cap[i]) for i=1:Nplot]
    y = [pa[2] + radius*sin(θ1cap[i]) for i=1:Nplot]
    plot!(plotObj,x,y)
end

# plots particles
# requires loading of Plots module
function plotParticle!(particle::Particle,plotObj::Plots.Plot,colorSym::Symbol,rm::Float64)
    θ = linspace(0,2*pi,40)
    xs = particle.xpos + 0.5*rm*cos.(θ)
    ys = particle.ypos + 0.5*rm*sin.(θ)

    plot!(plotObj,xs,ys)
end

end # SinkholePlotModule
