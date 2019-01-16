# Module used to Plot various aspects of
# Kept in seperate module to avoid loading time of Plots during normal
#   run-time

module partViz

# import only necessary functions from eFEMpart
using eFEMpart

using Plots, LaTeXStrings

#font = Plots.font("DejaVu Sans", 14)
#pyplot(size=(800,800),border=true,
#        guidefont=font, xtickfont=font, ytickfont=font, legendfont=font,
#        markersize=8,linewidth=3,ratio=:equal,grid=false)

# functions from partViz
export plotPoints,
       plotWall,
       plotWalls,
       plotParticles!
# functions from Plots
export @animate,
       plot,
       gif,
       Animation,
       frame,
       pyplot

function plotWalls(wallList::Vector{T},rm::Float64) where T<:AbstractWall
  N = length(wallList)

  p = plot(grid=false) #initialize plot
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
    θWALL = atan(y2-y1,x2-x1)

    Nθ = 15
    θ1 = range(-pi/2 + θWALL,stop=pi/2 + θWALL,length=Nθ)
    θ2 = range( pi/2 + θWALL,stop = 3*pi/2 + θWALL,length=Nθ)

    len = sqrt((x1-x2)^2 + (y1-y2)^2)
    
    tx = (x2-x1)/len
    ty = (y2-y1)/len

    nx = -ty
    ny = tx

    # 1 arc
    x = [x2 + radius*cos(θ1[i]) for i=1:Nθ]
    y = [y2 + radius*sin(θ1[i]) for i=1:Nθ]
    plot!(plotObj,x,y,label="",color=:black)

    # 2 line
    x = [x2 + radius*nx,x1 + radius*nx]
    y = [y2 + radius*ny,y1 + radius*ny]
    plot!(plotObj,x,y,label="",color=:black)

    # 3 arc 
    x = [x1 + radius*cos(θ2[i]) for i=1:Nθ]
    y = [y1 + radius*sin(θ2[i]) for i=1:Nθ]
    plot!(plotObj,x,y,label="",color=:black)
    
    # 4 line
    x = [x1 - radius*nx,x2 - radius*nx]
    y = [y1 - radius*ny,y2 - radius*ny]
    plot!(plotObj,x,y,label="",color=:black)
end
function plotWall(wall::CircleWall,plotObj::Plots.Plot,rm)
    radius = 0.5*rm
    Nplot = 50; Δθ = 2*pi/(Nplot-1)
    curvex = [wall.center.x + (wall.radius + radius)*cos((i-1)*Δθ) for i=1:Nplot]
    curvey = [wall.center.y + (wall.radius + radius)*sin((i-1)*Δθ) for i=1:Nplot]
    plot!(plotObj,curvex,curvey,label="",color=:black)
end
function plotWall(wall::ArcWall,plotObj::Plots.Plot,rm)
    radius = 0.5*rm

    pa = [wall.nodes[1].x,wall.nodes[1].y]
    c  = [wall.nodes[2].x,wall.nodes[2].y]
    pb = [wall.nodes[3].x,wall.nodes[3].y]
    r = sqrt((pa[1]-c[1])^2 + (pa[2]-c[2])^2)
    θ1 = atan(pa[2]-c[2],pa[1]-c[1])
    θ2 = atan(pb[2]-c[2],pb[1]-c[1])

    if θ1 > θ2
        θ2 += 2*pi
    end

    Nplot = 50; Δθ = (θ2 - θ1)/(Nplot-1)

    # 1 big arc -- done
    curvex = [c[1] + (r+radius)*cos((i-1)*Δθ+θ1) for i=1:Nplot]
    curvey = [c[2] + (r+radius)*sin((i-1)*Δθ+θ1) for i=1:Nplot]
    plot!(plotObj,curvex,curvey,label="",color=:black)

    # 2 far cap
    θWALL = θ2 + pi/2
    θ2cap = range(-pi/2 + θWALL, stop=pi/2 + θWALL,length=Nplot)
    x = [pb[1] + radius*cos(θ2cap[i]) for i=1:Nplot]
    y = [pb[2] + radius*sin(θ2cap[i]) for i=1:Nplot]
    plot!(plotObj,x,y,label="",color=:black)

    # 3 small arc 
    curvex = [c[1] + (r-radius)*cos((i-1)*Δθ+θ1) for i=1:Nplot]
    curvey = [c[2] + (r-radius)*sin((i-1)*Δθ+θ1) for i=1:Nplot]
    plot!(plotObj,curvex,curvey,label="",color=:black)
    
    # 4 close cap
    θWALL = θ1 + pi/2
    θ1cap = range( pi/2 + θWALL,stop = 3*pi/2 + θWALL,length=Nplot)
    x = [pa[1] + radius*cos(θ1cap[i]) for i=1:Nplot]
    y = [pa[2] + radius*sin(θ1cap[i]) for i=1:Nplot]
    plot!(plotObj,x,y,label="",color=:black)
end

# plots particles
# requires loading of Plots module
function plotParticle!(particle::Particle,plotObj::Plots.Plot,colorSym::Symbol,rm::Float64)
    θ = range(0,stop=2*pi,length=40)
    xs = particle.xpos .+ 0.5*rm*cos.(θ)
    ys = particle.ypos .+ 0.5*rm*sin.(θ)

    plot!(plotObj,xs,ys,label="",color=:black)
end

end # partViz
